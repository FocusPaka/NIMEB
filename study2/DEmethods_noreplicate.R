


#-----------------------------------------------------#
#---DE methods used in data without replicate, include 
#---MEB, HTN, edgeR, Lib Size, DESeq, NOISeq (6 methods)
#-----------------------------------------------------#

runDE <- 
  function(countsTable, condsAB, run_MEB=TRUE, train_id, 
           gamma=seq(0.0001,0.005,0.0001), nu=0.01, reject_rate=0.05,
           run_HTN=TRUE, run_edgeR=TRUE, run_LibSize=TRUE, 
           run_DESeq=TRUE, run_NOISeq=TRUE){
    
    #--------
    #--MEB
    #--------
    if(run_MEB){
      resMEB <- myMEB(countsTable, train_id, gamma, nu, reject_rate)
    }
  
    #--------
    #--HTN
    #--------
    if(run_HTN){
      resHTN <- myHTN(countsTable,train_id)
    }
    
    #--------
    #--edgeR
    #--------
    if(run_edgeR){
      resedgeR <- myedgeR(countsTable)
    }
    
    #-----------
    #--Lib Size
    #-----------
    if(run_LibSize){
      resLibSize <- myLibSize(countsTable)
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
    
    pfull <- data.frame(HTNseq=ifelse(rep(run_HTN, ng), resHTN, NA),
                         edgeR=ifelse(rep(run_edgeR, ng), resedgeR, NA),
                         LibSize=ifelse(rep(run_LibSize, ng), resLibSize, NA),
                         DESeq=ifelse(rep(run_DESeq, ng), resDESeq, NA),
                         NOISeq=ifelse(rep(run_NOISeq, ng), 1-resNOISeq, NA)
                         )
    
    #list(MEBmodel=resMEB, pfull=pfull)
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
                 nu = nu, tolerance = 0.001, 
                 shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                 na.action = na.omit)
    train_error[k] <- 1 - model$tot.accuracy/100
  }
  
  gamma_num_new <- which.min(abs(train_error-reject_rate))
  model_new <- svm(trainData, y = NULL, scale = FALSE, type = "one-classification", 
                   kernel = "radial", gamma = gamma[gamma_num_new],
                   nu = nu, tolerance = 0.001, 
                   shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                   na.action = na.omit)
  
  #return(model=model_new)
  list(model=model_new,gamma = gamma[gamma_num_new])
}


#---HTN method---#
myHTN <- function(countsTable,train_id){
  
  print("run HTN.")
  library(statmod)           #sage.test() funciton in this package
  source("other_functions.R")  
  
  hkeep <- train_id
  sc_factor <- SCBNM(countsTable,hkeep,a=0.05)
  pval <- sage.test(countsTable[,1], countsTable[,2], n1=sum(countsTable[,1])/sqrt(sc_factor), 
                     n2=sum(countsTable[,2])*sqrt(sc_factor))
  
  return(pval=pval)
}



#---edgeR method(TMM normalization)---#
myedgeR <- function(countsTable){
  
  print("run edgeR")
  #library(edgeR)
  
  fW <- calcNormFactors_new(countsTable,logratioTrim=0.3, sumTrim=0.05)[2]
  sM <- sage.test(countsTable[,1], countsTable[,2], n1=sum(countsTable[,1])/sqrt(fW), 
                  n2=sum(countsTable[,2])*sqrt(fW))
  
  return(pval=sM)
}



#---Lib Size method---#
myLibSize <- function(countsTable){
  
  print("run Lib Size.")
  
  s <- sage.test(countsTable[,1], countsTable[,2], n1=sum(countsTable[,1]), 
                 n2=sum(countsTable[,2]))
  
  return(pval=s)
}



#---DESeq method---#
#myDESeq <- function(countsTable,condsAB){
  
#  print("run DESeq.")
#  library(DESeq)
  
#  type <- factor(condsAB)
#  cds <- newCountDataSet(countsTable,type)
  
  #对于没有生物学重复
#  cds <- estimateSizeFactors(cds)
#  cds <- estimateDispersions(cds, method="blind", fitType="local", sharingMode="fit-only" )
#  res <- nbinomTest(cds,condsAB[1], condsAB[2])
#  pval <- res$pval
 
#  return(pval=pval)
#}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DESeq
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
myDESeq = function(countsTable, condsAB, colAB=NULL, method=NULL, 
                    sharingModeAB="fit-only"){
  
  library(DESeq)
  print("run DESeq.")  
  
  #filt=filter1(countsTable)
  #countsTable=filt$data
  #keep=filt$keep
  #countsTable[countsTable==0.5]=1
  
  
  factors <- unique(condsAB)
  
  if(is.null(colAB)){colAB=1:ncol(countsTable);}
  cds = newCountDataSet(countsTable, condsAB); 
  cds = estimateSizeFactors( cds );
  if(length(condsAB)>2){
    if(is.null(method)){
      cds = estimateDispersions( cds, fitType="local"); 
    } else {
      cds = estimateDispersions(cds, method=method, fitType="local", sharingMode = "maximum"); 
    }
  } else {
    #when no replicates.
    cds = estimateDispersions(cds, method="blind",fitType="local", 
                              sharingMode= sharingModeAB);    
  }
  res = nbinomTest( cds, factors[[1]], factors[[2]] );  
  pval <- res$pval
  
  return(pval=pval)  
}





#---NOISeq method---#
myNOISeq <- function(countsTable,condsAB){
  
  print("run NOISeq.")
  library(NOISeq)
  
  myfactors = data.frame(condsAB = condsAB)
  myData=NOISeq::readData(data = countsTable, factors = myfactors)
  mynoiseq=noiseq(myData,k=0.5,norm="tmm",factor="condsAB",pnr=0.2,nss=5,v=0.02,
                  lc=1,replicates="no")
  
  noiseq_prob <- mynoiseq@results[[1]][,"prob"]
  
  return(pval=noiseq_prob)
}











