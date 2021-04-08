

#-----------------------------------------------------------#
#---DE methods used in data with replicates, include 
#---MEB, Marioni(no norm), Marioni (norm), edgeR(no norm),
#---edgeR(norm), cloonan, HTN, DESeq, NOISeq(9 methods)
#-----------------------------------------------------------#

runDEs <- 
  function(countsTable, condsAB, run_MEB=TRUE, train_id, 
           gamma=seq(0.0001,0.005,0.0001), nu=0.01, reject_rate=0.05,
           run_Marioni0=TRUE, run_Marioni=TRUE,genelength, run_edgeR0=TRUE, 
           run_edgeR=TRUE,run_cloonan=TRUE,run_HTN=TRUE, 
           run_DESeq=TRUE, run_NOISeq=TRUE){
    
    #--------
    #--MEB
    #--------
    if(run_MEB){
      resMEB <- myMEB(countsTable, train_id, gamma, nu, reject_rate)
    }
    
    #-----------------------
    #--Marioni(no normalize/normalize)
    #-----------------------
    if(run_Marioni){
      resMarioni <- myMarioni(countsTable,condsAB)
    }
    
    
    #-----------------------
    #--edgeR(no normalize/normalize)
    #-----------------------
    if(run_edgeR){
      resedgeR <- myedgeR(countsTable,condsAB)
    }
  
    #-----------
    #--cloonan
    #-----------
    if(run_cloonan){
      rescloonan <- mycloonan(countsTable,condsAB, genelength)   
    }
    
    
    #--------
    #--HTN
    #--------
    if(run_HTN){
      resHTN <- myHTN(countsTable,condsAB,train_id)
    }
    
    
    #--------
    #--DESeq
    #--------
    if(run_DESeq){
      resDESeq <- myDESeq(countsTable,condsAB)
    }
    

    #--------
    #--NOISeq
    #--------
    if(run_NOISeq){
      resNOISeq <- myNOISeq(countsTable,condsAB)
    }
    
    #-------------------
    # Combine the "p-values" from the seven methods into a data frame
    #--------------------
    
    ng <- nrow(countsTable)
    pfull <- data.frame(Marioni0=ifelse(rep(run_Marioni0, ng), resMarioni$pval0, NA),
                        Marioni=ifelse(rep(run_Marioni, ng), resMarioni$pval, NA),
                        edgeR0=ifelse(rep(run_edgeR0, ng), resedgeR$pval0, NA),
                        edgeR=ifelse(rep(run_edgeR, ng), resedgeR$pval, NA),
                        cloonan=ifelse(rep(run_cloonan, ng), rescloonan, NA),
                        HTNseq=ifelse(rep(run_HTN, ng), resHTN, NA),
                        DESeq=ifelse(rep(run_DESeq, ng), resDESeq, NA),
                        NOISeq=ifelse(rep(run_NOISeq, ng), 1-resNOISeq, NA)
    )
    
    list(MEBmodel=resMEB$model,gamma = resMEB$gamma, pfull=pfull)
  }



#---MEB method---#
myMEB <- function(countsTable, train_id, gamma, nu, reject_rate){
  
  print("run MEB.")
  library("e1071")
  
  trainData <- countsTable[train_id,]
  train_error <- numeric(length(gamma))
  for(k in 1:length(gamma)){
    model <- svm(trainData, y = NULL, scale = FALSE, type = "one-classification", 
                 kernel = "radial", gamma = gamma[k],
                 nu = 0.01, tolerance = 0.001, 
                 shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                 na.action = na.omit)
    train_error[k] <- 1 - model$tot.accuracy/100
  }
  
  gamma_num_new <- which.min(abs(train_error-reject_rate))
  model_new <- svm(trainData, y = NULL, scale = FALSE, type = "one-classification", 
                   kernel = "radial", gamma = gamma[gamma_num_new],
                   nu = 0.01, tolerance = 0.001, 
                   shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                   na.action = na.omit)
  
  list(model=model_new,gamma = gamma[gamma_num_new])
}


#----Marioni(no normalize/normalize)----#
myMarioni <- function(countsTable,condsAB){
  
  print("run Marioni.")
  gg <- condsAB
  gl <- levels(condsAB)
  pm <- Poisson.model.new(countsTable, which(gg==gl[1]), which(gg==gl[2]), calcFactor=FALSE)
  pmAdj <- Poisson.model.new(countsTable, which(gg==gl[1]), which(gg==gl[2]), calcFactor=TRUE)
  
  list(pval0=pm$stats$pval,pval=pmAdj$stats$pval) 
}


#-----edgeR(no normalize/normalize)---#
myedgeR <- function(countsTable,condsAB){
  
  print("run edgeR.")
  gg <- condsAB
  gl <- levels(condsAB)
  fW <- calcNormFactors_new(countsTable)
  lambda <- rowSums(countsTable)/sum(countsTable)
  expMean <- outer(lambda, colSums(countsTable))
  cs <- colSums(countsTable)
  effM <- cs*fW
  expMeanAdj <- outer(lambda, effM)
  
  exactP <- exactTestPoisson(countsTable, group1Ind=which(gg==gl[1]), 
                             group2Ind= which(gg==gl[2]), 
                             expMean, verbose = TRUE)
  exactPadj <-exactTestPoisson(countsTable, group1Ind=which(gg==gl[1]), 
                               group2Ind=which(gg==gl[2]), 
                               expMeanAdj, verbose = TRUE)
  
  list(pval0=exactP,pval=exactPadj)
}


#myedgeR = function(countsTable, condsAB, colAB=NULL, sfact=NULL){

#  library(edgeR)
#  print("run edgeR.")

#filt=filter2(countsTable)
#countsTable=filt$data
#keep=filt$keep

#  if(is.null(colAB)){colAB=1:ncol(countsTable);}      
#  dgl = DGEList( counts=countsTable[,colAB], group=condsAB)
#  if(is.null(sfact)){
#    dgl = calcNormFactors(dgl);
#  } else {
#    dgl$samples$norm.factors = sfact;   
#  }
# print(dgl$common.dispersion);
# print(sfact);
#  dgl = estimateCommonDisp(dgl);  
#  dgl = estimateTagwiseDisp(dgl); 
#  edgerRes.tw = exactTest( dgl , dispersion="auto");  
#padj= p.adjust(edgerRes.tw$table$PValue, method="BH")
#print(dgl$samples);
#  pval=padj=rep(1,nrow(countsTable))
#  pval=edgerRes.tw$table$PValue
#  padj <- p.adjust(edgerRes.tw$table$PValue, method="BH")
#  list(pval = pval, padj=padj);      
#}







#------cloonan-------#
mycloonan <- function(countsTable,condsAB, genelength){
  
  print("run cloonan.")
  d <- log2(countsTable+1)/outer(genelength,rep(1,ncol(countsTable)))
  d <- normalizeQuantiles(d)
  
  mm <- model.matrix(~condsAB)   
  f <- lmFit(d, mm)
  f <- eBayes(f)
  cloonanpval <- f$p.value[,2]
  
  return(pval=cloonanpval)
}



#------HTN method----#
myHTN <- function(countsTable,condsAB,train_id){
  
  print("run HTN.")
  gg <- condsAB
  gl <- levels(condsAB)
  hkeep <- train_id
  nn<-length(countsTable[1,])
  fww<-rep(1,nn)
  fww[1] = 1
  for(t in 2:nn){
    fww[t]<-SCBNM(countsTable[,c(1,t)],hkeep,a=0.05)
  }
  lambda2 <- rowSums(countsTable)/sum(countsTable)
  cs <- colSums(countsTable)
  effM3 <- cs*fww
  expMeanAdj3 <- outer(lambda2, effM3)
  exactPadj3 <-exactTestPoisson(countsTable, group1Ind=which(gg==gl[1]),
                                group2Ind= which(gg==gl[2]),
                                expMeanAdj3, verbose = TRUE)
  return(pval=exactPadj3)
}



#-----DESeq method----#
myDESeq <- function(countsTable,condsAB){
  
  print("run DESeq.")
  library("DESeq2")
  exprSet=countsTable
  group_list=condsAB
  
  colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds2 <- DESeq(dds)  
  res <-  results(dds2, contrast=c("group_list","1","2"))
  pval=res$pvalue
  
  return(pval=pval)
}



#-----NOISeq method----#
myNOISeq <- function(countsTable,condsAB){
  
  print("run NOISeq.")
  library("NOISeq")
  
  myfactors = data.frame(condsAB = condsAB)
  myData=NOISeq::readData(data = countsTable, factors = myfactors)
  mynoiseq=noiseq(myData,k=0.5,norm="tmm",factor="condsAB",pnr=0.2,nss=5,v=0.02,
                  lc=1,replicates="biological")
  
  noiseq_prob <- mynoiseq@results[[1]][,"prob"]
  
  return(pval=noiseq_prob)
}















