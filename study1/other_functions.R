#(1) SBNM method to choose normalization factor.
SCBNM<-function(xx,housekeep,a=0.05)
{
  library(edgeR)
  library(statmod)
  fW <- calcNormFactors_new(xx,logratioTrim=0.3, sumTrim=0.05)[2]
  fW4<-seq(max(fW-0.5,0.01),fW+0.5,0.1)
  
  n<-length(fW4)
  fdr<-rep(0,n)
  for(j in 1:n){
    sMm4<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW4[j]), 
                     n2=sum(xx[,2])*sqrt(fW4[j]))
    fdr[j]<-sum(sMm4[housekeep]<a)
  }
  fw5<-fW4[which.min(abs(fdr-a))]
  fW41<-seq(max(fw5-0.25,0.01),fw5+0.25,0.005)
  n1<-length(fW41)
  fdr1<-rep(0,n1)
  for(j in 1:n1){
    sMm41<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW41[j]), 
                      n2=sum(xx[,2])*sqrt(fW41[j]))
    fdr1[j]<-sum(sMm41[housekeep]<a)
    if(j %% 10==0) cat(".")
  }
  fw51<-fW41[which.min(abs(fdr1-a))]
  return(fw51)
}


calcNormFactors_new <- function(dataMatrix, refColumn=1, logratioTrim=.3, 
                                sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  if( !is.matrix(dataMatrix) )
    stop("'dataMatrix' needs to be a matrix")
  if( refColumn > ncol(dataMatrix) )
    stop("Invalid 'refColumn' argument")
  apply(dataMatrix,2,.calcFactorWeighted,ref=dataMatrix[,refColumn], logratioTrim=logratioTrim, 
        sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
}

.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, 
                                sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  
  if( all(obs==ref) )
    return(1)
  
  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}








sage.test2 <- function(x, y, n1=sum(x), n2=sum(y))
{
  if(any(is.na(x)) || any(is.na(y))) stop("missing values not allowed")
  nn<-sum(x)*sqrt(n1/n2)
  mm<-sum(y)*sqrt(n2/n1)
  x <- round(x)
  y <- round(y)
  if(any(x<0) || any(y<0)) stop("x and y must be non-negative")
  if(length(x) != length(y)) stop("x and y must have same length")
  n1 <- round(n1)
  n2 <- round(n2)
  #if(!missing(n1) && any(x>n1)) stop("x cannot be greater than n1")
  #if(!missing(n2) && any(y>n2)) stop("y cannot be greater than n2")
  size <- x+y
  p.value <- rep(1,length(x))
  if(n1==n2) {
    i <- (size>0)
    if(any(i)) {
      x <- pmin(x[i],y[i])
      size <- size[i]
      p.value[i] <- pbinom(x,size=size,prob=0.5)+pbinom(size-x+0.5,size=size,
                                                        prob=0.5,lower.tail=FALSE)
    }
    return(p.value)
  }
  prob <- n1/(n1+n2)
  nn <- round(nn)
  mm <- round(mm)
  
  if(any(big <- size>10000)) {
    ibig <- (1:length(x))[big]
    for (i in ibig) p.value[i] <- chisq.test(matrix(c(x[i],y[i],(nn-x[i]),
                                                      (mm-y[i])),2,2))$p.value
  }
  size0 <- size[size>0 & !big]
  if(length(size0)) for (isize in unique(size0)) {
    i <- (size==isize)
    p <- dbinom(0:isize,p=prob,size=isize)
    o <- order(p)
    cumsump <- cumsum(p[o])[order(o)]
    p.value[i] <- cumsump[x[i]+1]
  }
  p.value
}





Poisson.model.new <- function(countMatrix,group1,group2, ref=1, calcFactor=TRUE){
  
  Poisson.glm.pval <- vector()
  Fold.changes <- vector()
  
  props <- countMatrix / outer( rep(1,nrow(countMatrix)), colSums(countMatrix) )
  
  refS <- colSums(countMatrix[,c(group1,group2)])
  
  if( calcFactor ) {
    require(edgeR)
    CS <- calcNormFactors(countMatrix[,c(group1,group2)])
  } else {
    CS <- rep(1,length(group1)+length(group2))
  }
  
  offsets <- log(CS)+log(refS)
  
  sample.f <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
  
  for (i in 1:(nrow(countMatrix))){
    S1 <- countMatrix[i,group1] 
    S2 <- countMatrix[i,group2] 
    In <- c(S1,S2)
    In <- as.vector(unlist(In))
    GLM.Poisson <- glm(In ~ 1 + sample.f + offset(offsets),family=poisson)
    Poisson.glm.pval[i] <- anova(GLM.Poisson,test="Chisq")[5][2,1]
    Fold.changes[i] <- exp(GLM.Poisson$coefficients[1])/(exp(GLM.Poisson$coefficients[1]+GLM.Poisson$coefficients[2]))
    if(i %% 100==0) cat(".")
  }
  cat("\n")
  
  #output <- matrix(ncol=2,nrow=nrow(countMatrix))
  #output[,1] <- Poisson.glm.pval
  #output[,2] <- Fold.changes
  #output <- as.data.frame(output)
  #names(output) <- c("pval","FC")
  
  list(stats=data.frame(pval=Poisson.glm.pval, FC=Fold.changes),offsets=offsets,factors=CS)
}



exactTestPoisson <- function(dataMatrix, meanMatrix, group1Ind, group2Ind, verbose=TRUE) {
  if(length(group1Ind) < 2){
    y1 <- dataMatrix[,group1Ind]
    y2 <- dataMatrix[,group2Ind]
    m1 <- meanMatrix[,group1Ind]
    m2 <- meanMatrix[,group2Ind]
  }else{
    y1 <- rowSums(dataMatrix[,group1Ind])
    y2 <- rowSums(dataMatrix[,group2Ind])
    m1 <- rowSums(meanMatrix[,group1Ind])
    m2 <- rowSums(meanMatrix[,group2Ind])
  }
  
  N <- rowSums( dataMatrix[,c(group1Ind,group2Ind)] )
  
  pvals <- rep(NA, nrow(dataMatrix))
  
  for (i in 1:length(pvals)) {
    v <- 0:N[i]
    p.top <- dpois(v, lambda=m1[i]) * dpois(N[i]-v, lambda=m2[i])
    p.obs <- dpois(y1[i], lambda=m1[i]) * dpois(y2[i], lambda=m2[i])
    p.bot <- dpois(N[i], lambda=m1[i]+m2[i])
    keep <- p.top <= p.obs
    pvals[i] <- sum(p.top[keep]/p.bot)
    if(N[i]>10000) {	
      pvals[i]<- chisq.test(matrix(c(y1[i],y2[i],(sum(y1)-y1[i]),(sum(y2)-y2[i])),2,2))$p.value
    }
    if (verbose)
      if (i%%1000 == 0)
        cat(".")
  }
  if (verbose)
    cat("\n")  
  pvals
  
}




##############################################
kernel <- function(x,sv,gamma){
  k <- numeric(nrow(sv))
  for(i in 1:nrow(sv)){
    k[i] <- exp(-gamma*sum((sv[i,]-x)^2))
  }
  return(k)
}



##############################################
decision_function <- function(x,model,gamma){
  index <- model$coefs
  sv <- model$SV
  dec_val <- sum(index*kernel(x,sv=sv,gamma))
  return(dec_val)
}























