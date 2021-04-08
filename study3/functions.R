#生成数据函数
#没有重复样本的情况
generateDataset <- function(commonTags=15000, uniqueTags=c(1000,3000), 
                            group=c(1,2), libLimits=c(.9,1.1)*1e6, 
                            empiricalDist=NULL, randomRate=1/100, 
                            pDifferential=.05, pUp=.5, foldDifference=2) {
    
    # some checks
    group <- as.factor(group)
    stopifnot( length(group) == length(uniqueTags) )
    stopifnot( nlevels(group) == 2 ) 
    
    # define where to take random sample from (empirical distribution OR random exponential)
    if(is.null(empiricalDist))
        exampleCounts <- ceiling(rexp(commonTags,rate=randomRate))
    else
        exampleCounts <- empiricalDist
    
    exampleLambda <- exampleCounts/sum(exampleCounts)
    
    # set up libraries
    nLibraries <- length(uniqueTags)
    libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )
    
    # vector of starts/stops for the unique Tags
    en <- commonTags + cumsum(uniqueTags)
    st <- c(commonTags+1,en[-nLibraries]+1)
    
    # create matrix of LAMBDA(=relative expression levels)
    LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
    LAMBDA[1:commonTags,] <- sample(exampleLambda, commonTags, replace=TRUE)
    
    # set unique tag totals
    for(i in 1:nLibraries)
        if(uniqueTags[i] > 0)
            LAMBDA[st[i]:en[i],i] <- sample(exampleLambda, uniqueTags[i])    
    
    ind <- seq_len(floor(pDifferential*commonTags))
    g <- group == levels(group)[1]
    if(length(ind)>0) {
        fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
        LAMBDA[ind,g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*fcDir)
        LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
    }
    
    sampFactors <- colSums(LAMBDA)
    
    sampFactorsM <- outer(rep(1,max(en)),sampFactors)
    libSizesM <- outer(rep(1,max(en)),libSizes)
    
    # create observed means
    MEAN <- LAMBDA / sampFactorsM * libSizesM
    
    true.log2foldchange <- log2(LAMBDA[1:commonTags,2]/LAMBDA[1:commonTags,1])
    
    # sample observed data (column sums will be *close* to set library sizes)
    DATA <- matrix(rpois(length(MEAN), lambda=MEAN),ncol=nLibraries)
    
    #trueFactors <- colSums(MEAN[1:commonTags,])
    trueFactors <- colSums(MEAN[(length(ind)+1):commonTags,])
    trueFactors <- trueFactors/trueFactors[1]
    
    list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, 
         trueFactors=trueFactors, 
         truelog2foldchanges=true.log2foldchange,
         group=group, libSizes=libSizes,  
         differentialInd=c(ind,(commonTags+1):nrow(DATA)),
         commonInd=1:commonTags)
}




##################################
#有重复样本的情况
generateDataset2 <- function(commonTags=15000, uniqueTags=c(1000,3000), 
                             group=c(1,2), libLimits=c(.9,1.1)*1e6, 
                             empiricalDist=NULL, lengthDist=NULL, 
                             pDifferential=.05, pUp=.5, foldDifference=2, 
                             nreps=c(2,2)) {
    
    # some checks
    stopifnot( length(group) == length(uniqueTags) )
    stopifnot( length(group) == length(nreps) )
    stopifnot( length(empiricalDist) == length(lengthDist) )
    group <- as.factor(rep(group,nreps))
    stopifnot( nlevels(group) == 2 ) 
    
    print(group)
    
    #exampleCounts <- empiricalDist/lengthDist
    exampleCounts <- empiricalDist
    exampleLambda <- exampleCounts/sum(exampleCounts)
    exampleIds <- seq_len( length(empiricalDist) )
    
    # set up libraries
    nLibraries <- sum( nreps )
    libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )
    
    # vector of starts/stops for the unique Tags
    en <- commonTags + cumsum(uniqueTags)
    st <- c(commonTags+1,en[-nLibraries]+1)
    
    # create matrix of LAMBDA(=relative expression levels)
    LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
    
    ID <- rep(0, max(en))
    ID[1:commonTags] <- sample(exampleIds, commonTags, replace=TRUE)
    LAMBDA[1:commonTags,] <- exampleLambda[ ID[1:commonTags] ]
    
    # set unique tag totals
    for(i in 1:length(nreps))
        if(uniqueTags[i] > 0) {
            ID[st[i]:en[i]] <- sample(exampleIds, uniqueTags[i], replace=TRUE)
            LAMBDA[st[i]:en[i],group==levels(group)[i]] <- exampleLambda[ ID[st[i]:en[i]] ]
        }
    
    
    # set fold tags 
    g <- group == levels(group)[1]
    ind <- seq_len(floor(pDifferential*commonTags))
    if(length(ind)>0) {
        fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
        LAMBDA[ind,g] <- LAMBDA[ind,g]*exp(log(foldDifference)/2*fcDir)
        LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
    }
    
    sampFactors <- colSums(LAMBDA)
    
    sampFactorsM <- outer( rep(1,max(en)), sampFactors )    #rep(1,max(en))行，sampFactors列的矩阵
    libSizesM <- outer(  rep(1,max(en)), libSizes )         #rep(1,max(en))行，libSizes列的矩阵
    
    # create observed means
    MEAN <- LAMBDA / sampFactorsM * libSizesM  # to get the totals to sum to 1
    
    # sample observed data (column sums will be *close* to set library sizes)
    DATA <- matrix(0, nrow=nrow(LAMBDA), ncol=nLibraries)
    DATA <- matrix(rpois(length(MEAN), lambda=MEAN), ncol=nLibraries)    #按列求
    
    trueFactors <- colSums(MEAN[1:commonTags,])   
    trueFactors <- trueFactors/trueFactors[1]
    
    true.log2foldchange <- log2(LAMBDA[1:commonTags,3]/LAMBDA[1:commonTags,1])
    
    colnames(DATA) <- paste(paste("group",group,sep=""),1:ncol(DATA),sep=".")
    
    list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, trueFactors=trueFactors, truelog2foldchanges = true.log2foldchange,group=group, libSizes=libSizes,  
         differentialInd=c(ind,(commonTags+1):nrow(DATA)), commonInd=1:commonTags, ID=ID, length=lengthDist[ID])
}


#########################################
takeSubset <- function(obj, subsetInd) {
    allInd <- 1:nrow(obj$DATA)
    commonInd <- allInd %in% obj$commonInd
    differentialInd <- allInd %in% obj$differentialInd
    list(DATA=obj$DATA[subsetInd,], LAMBDA=obj$LAMBDA[subsetInd,], 
         #trueFactors=obj$trueFactors, 
         group=obj$group, 
         libSizes=obj$libSizes, 
         differentialInd=which(differentialInd[subsetInd]), 
         commonInd=which(commonInd[subsetInd]),
         ID=obj$ID[subsetInd], length=obj$length[subsetInd])
}





################################################
select_conser4 <- function(data, comm, dii, prop = c(1,2), n = 150){
    extreme <- apply(data[comm,], 1, function(x) max(x) - min(x))
    select <-  (extreme < n) & (extreme > 5)
    count1 <- which(select == TRUE)
    
    temp <- (comm - length(dii))%in% count1
    #ratio <- numeric(length(comm))
    m <- ncol(data)
    ratio <- apply(data[comm,], 1, function(x) max(mean(x[1:(m/2)]),
                                                   mean(x[(m/2+1):m]))/min(mean(x[1:(m/2)]),
                                                                           mean(x[(m/2+1):m])))
    
    count2 <- intersect(which((ratio >= 1) & (ratio <= 2)),which(!temp))
    count <- union(count1,count2)
    
    select_com <- count + length(dii)
    return(select_com)
}








##########################################################
#不同方法类函数
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



##################################
#(1) SBNM method to choose normalization factor.
SCBNM<-function(xx,housekeep,a=0.05)
{
    library(edgeR)
    library(statmod)
    fW <- calcNormFactors_new(xx,logratioTrim=0.3, sumTrim=0.05)[2]
    fW4<-seq(max(0.25,fW-0.5),fW+0.5,0.1)
    
    n<-length(fW4)
    fdr<-rep(0,n)
    for(j in 1:n){
        sMm4<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW4[j]), 
                         n2=sum(xx[,2])*sqrt(fW4[j]))
        fdr[j]<-sum(sMm4[housekeep]<a)/length(housekeep)
    }
    fw5<-fW4[which.min(abs(fdr-a))]
    fW41<-seq(max(0.25,fw5-0.25),fw5+0.25,0.005)
    n1<-length(fW41)
    fdr1<-rep(0,n1)
    for(j in 1:n1){
        sMm41<-sage.test2(xx[,1], xx[,2], n1=sum(xx[,1])/sqrt(fW41[j]), 
                          n2=sum(xx[,2])*sqrt(fW41[j]))
        fdr1[j]<-sum(sMm41[housekeep]<a)/length(housekeep)
        if(j %% 10==0) cat(".")
    }
    fw51<-fW41[which.min(abs(fdr1-a))]
    return(fw51)
}

#####################################
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





###################################
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


##################################
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


##################################
exactTestPoisson2 <- function(dataMatrix, meanMatrix, group1Ind, group2Ind, verbose=TRUE) {
    
    y1 <- rowSums(dataMatrix[,group1Ind])
    y2 <- rowSums(dataMatrix[,group2Ind])
    m1 <- rowSums(meanMatrix[,group1Ind])
    m2 <- rowSums(meanMatrix[,group2Ind])
    
    N <- rowSums( dataMatrix[,c(group1Ind,group2Ind)] )
    
    pvals <- rep(NA, nrow(dataMatrix))
    ind <- rep(0, nrow(dataMatrix))
    N1<-sum(dataMatrix[,group1Ind])
    N2<-sum(dataMatrix[,group2Ind])
    for (i in 1:length(pvals)) {
        v <- 0:N[i]
        p.top <- dpois(v, lambda=m1[i]) * dpois(N[i]-v, lambda=m2[i])
        p.obs <- dpois(y1[i], lambda=m1[i]) * dpois(y2[i], lambda=m2[i])
        p.bot <- dpois(N[i], lambda=m1[i]+m2[i])
        if(ppois(y1[i], lambda=m1[i])<ppois(y2[i], lambda=m2[i])) {ind[i]<-1}
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
    list(pvals=pvals,ind=ind)
}




#####################################################
#训练，选择参数
train_procedure <- function(x_train, gamma, nu){
    train_error <- matrix(0, length(gamma),length(nu))
    train_nSV <- matrix(0,length(gamma),length(nu))
    for(i in 1:length(gamma)){
        for(j in 1:length(nu)){
            model <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification", 
                         kernel = "radial", gamma = gamma[i],
                         nu = nu[j], tolerance = 0.001, 
                         shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                         na.action = na.omit)
            
            train_error[i,j] <- 1 - model$tot.accuracy/100
            train_nSV[i,j] <- model$tot.nSV
        }
    }
    list(error=train_error,nSV=train_nSV)
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

