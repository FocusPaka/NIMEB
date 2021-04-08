

#-----------------------------------------------------------#
#---DE methods used in data with replicates, include 
#---NIMEB, Library size, edgeR, HTN, DESeq2, NOISeq (6 methods)
#-----------------------------------------------------------#

runDEs <- 
    function(countsTable, condsAB, run_NIMEB=TRUE, train_id, 
             gamma=seq(0.0001,0.005,0.0001), nu=0.01, reject_rate=0.05,
             run_librarySize=TRUE, run_edgeR=TRUE, run_HTN=TRUE, 
             run_DESeq2=TRUE, run_NOISeq=TRUE){
        
        #--------
        #--NIMEB
        #--------
        if (run_NIMEB) {
            resNIMEB <- myNIMEB(countsTable, train_id, gamma, nu, reject_rate)
        }
        
        
        #-----------------------
        #--Library size
        #-----------------------
        if (run_librarySize) {
            reslibrarySize <- mylibrarySize(countsTable,condsAB)
        }
        
        
        #-----------------------
        #--edgeR
        #-----------------------
        if (run_edgeR) {
            resedgeR <- myedgeR(countsTable,condsAB)
        }
        

        #--------
        #--HTN
        #--------
        if (run_HTN) {
            resHTN <- myHTN(countsTable,condsAB,train_id)
        }
        
        
        #--------
        #--DESeq2
        #--------
        if (run_DESeq2) {
            resDESeq <- myDESeq2(countsTable,condsAB)
        }
        
        
        #--------
        #--NOISeq
        #--------
        if (run_NOISeq) {
            resNOISeq <- myNOISeq(countsTable,condsAB)
        }
        
        
        #-------------------
        # Combine the "p-values" from the five methods into a data frame
        #--------------------
        ng <- nrow(countsTable)
        
        pfull <- data.frame(
            librarySize = ifelse(rep(run_librarySize, ng), 
                                 reslibrarySize, NA),
            edgeR = ifelse(rep(run_edgeR, ng), resedgeR$pval, NA),
            HTNseq = ifelse(rep(run_HTN, ng), resHTN, NA),
            DESeq = ifelse(rep(run_DESeq2, ng), resDESeq, NA),
            NOISeq = ifelse(rep(run_NOISeq, ng), 1 - resNOISeq, NA)
        )
        
        list(model = resNIMEB$model, pfull = pfull)
    }



#---NIMEB method---#
myNIMEB <- function(countsTable, train_id, gamma, nu, reject_rate){
    
    print("run NIMEB.")
    library("e1071")
    
    trainData <- countsTable[train_id,]
    train_error <- numeric(length(gamma))
    for (k in 1:length(gamma)) {
        model <- svm(trainData, y = NULL, scale = FALSE, type = "one-classification", 
                     kernel = "radial", gamma = gamma[k],
                     nu = 0.01, tolerance = 0.001, 
                     shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                     na.action = na.omit)
        train_error[k] <- 1 - model$tot.accuracy/100
    }
    
    gamma_num_new <- which.min(abs(train_error - reject_rate))
    
    model_new <- svm(trainData, y = NULL, scale = FALSE, type = "one-classification", 
                     kernel = "radial", gamma = gamma[gamma_num_new],
                     nu = 0.01, tolerance = 0.001, 
                     shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                     na.action = na.omit)
    
    list(model = model_new)
}


#----Library size----#
mylibrarySize <- function(countsTable,condsAB){
    
    print("run library Size.")
    gg <- condsAB
    libsize.pvalues <- Poisson.model(countMatrix = countsTable, 
                                     group1 = which(gg == unique(condsAB)[1]), 
                                     group2 = which(gg == unique(condsAB)[2]), 
                                     calcFactor = FALSE)
    
    return(pval = libsize.pvalues$stats$pval) 
}


#--------edgeR-------#

myedgeR = function(countsTable, condsAB){
    
    print("run edgeR.")
    library(edgeR)
    
    edgeR.dgelist <- DGEList(counts = countsTable, group = factor(condsAB))
    edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = "TMM")
    edgeR.dgelist <- estimateCommonDisp(edgeR.dgelist)
    edgeR.dgelist <- estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
    edgeR.test <- exactTest(edgeR.dgelist)
    edgeR.pvalues <- edgeR.test$table$PValue
    edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = "BH")
    
    list(pval = edgeR.pvalues, padj = edgeR.adjpvalues)      
}



#------HTN method----#
myHTN <- function(countsTable,condsAB,train_id){
    
    print("run HTN.")
    gg <- condsAB

    hkeep <- train_id
    nn <- ncol(countsTable)
    fww <- rep(1,nn)
    fww[1] <- 1
    for (t in 2:nn) {
        fww[t] <- SCBNM(countsTable[,c(1,t)], hkeep, a = 0.05)
    }
    
    lambda2 <- rowSums(countsTable)/sum(countsTable)
    cs <- colSums(countsTable)
    effM3 <- cs*fww
    meanMatrix <- outer(lambda2, effM3)
    htn.pvalues <-exactTestPoisson(countsTable, 
                                  group1Ind = which(gg == unique(condsAB)[1]),
                                  group2Ind = which(gg == unique(condsAB)[2]),
                                  meanMatrix = meanMatrix, verbose = TRUE)
    return(pval = htn.pvalues)
}



#-----------DESeq2----------------------------------
myDESeq2 <- function(countsTable, condsAB){
    
    library(DESeq2)
    print("run DESeq.")  
    
    conditions <- factor(condsAB)
    dds <- DESeqDataSetFromMatrix(countsTable, DataFrame(conditions),
                                  ~ conditions)

    dds <- DESeq(dds)
    res <- results(dds)
    
    return(pval = res$pvalue)   
}






#-----NOISeq method----#
myNOISeq <- function(countsTable, condsAB){
    print("run NOISeq.")
    library("NOISeq")
    rownames(countsTable) <- NULL
    
    myfactors <- data.frame(condsAB = condsAB)
    myData <- NOISeq::readData(data = countsTable, factors = myfactors)
    mynoiseq <- noiseq(myData,k=0.5,norm="tmm",factor="condsAB",pnr=0.2,
                       nss=5,v=0.02,lc=1,replicates="technical")
    
    noiseq_prob <- mynoiseq@results[[1]][,"prob"]

    return(prob = noiseq_prob)
}












