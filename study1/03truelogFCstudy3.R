rm(list=ls())


library(pROC)

source("functions.R")
#source("DEmethods_replicates.R")

#产生数据
data_primal <- read.table("Grimmond_lengths.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------
pD <- 0.3
nreps <- c(2,2)
xx <- generateDataset2(commonTags = 15000, uniqueTags=c(1000,800),
                       group=c(1,2), libLimits=c(.9,1.2)*1e6,
                       empiricalDist = data_primal$EB,
                       lengthDist = data_primal$transcript_length,
                       pDifferential = pD, pUp=0.9,
                       foldDifference= 4, nreps= nreps)



summary(xx$truelog2foldchanges)


truelogfoldchangesTotal <- xx$truelog2foldchanges

truelogfoldchangesTotalM <- matrix(0, nrow = length(truelogfoldchangesTotal),ncol = 2)

truelogfoldchangesTotalM[,1] <- truelogfoldchangesTotal
truelogfoldchangesTotalM[truelogfoldchangesTotal !=0,2] <- 1   #DE gene标志为1，non-DE标志为0

x <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==0)]    #先画non-DE
y <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==0),1]

x1 <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==1)]
y1 <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==1),1]
plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("true"," ",log[2],"FC")),main = "Study 3")
points(x1,y1,col="red",pch = 16,ylim = c(-3,3),xlim = c(0,15000))


#expression(paste(log[2],"FC=1")

















