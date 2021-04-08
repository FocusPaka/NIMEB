#多个replicates=2,3,4,5,6,7的情况（循环）
rm(list=ls())

#setwd("E:/paper/paper3_SVM/svm/code0510")
source("functions.R")
#source("DEmethods_replicates.R")

#产生数据
data_primal <- read.table("E:\\paper\\paper3_SVM\\svm\\code0510\\data\\Grimmond_lengths.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)





nreps <- c(2,2)
xx_new1 <- generateDataset2(commonTags=15000, uniqueTags=c(1000,800), 
                            empiricalDist = data_primal$EB,
                            lengthDist = data_primal$transcript_length,
                            pDifferential=0.6, foldDifference=4, pUp=0.9,
                            libLimits=c(.9,1.2)*1e6,nreps = nreps)

xx_new2 <- generateDataset2(commonTags = 10000, uniqueTags=c(2000,1000), 
                            libLimits=c(.9,1.2)*1e6,empiricalDist = data_primal$EB,
                            lengthDist = data_primal$transcript_length, 
                            pDifferential = 0.4, pUp=.1, 
                            foldDifference= 8,nreps = nreps)




# View(xx_new1$LAMBDA)
# 
# 
# hist(xx_new1$truelog2foldchanges)
# summary(xx_new1$truelog2foldchanges)
# 
# plot(xx_new1$truelog2foldchanges)
# 
# 
# hist(xx_new2$truelog2foldchanges)
# summary(xx_new2$truelog2foldchanges)
# 
# plot(xx_new2$truelog2foldchanges)




truelogfoldchangesTotal <- c(xx_new1$truelog2foldchanges,xx_new2$truelog2foldchanges)
summary(truelogfoldchangesTotal)
hist(truelogfoldchangesTotal)
plot(truelogfoldchangesTotal)



truelogfoldchangesTotalM <- matrix(0, nrow = length(truelogfoldchangesTotal),ncol = 2)

truelogfoldchangesTotalM[,1] <- truelogfoldchangesTotal
truelogfoldchangesTotalM[truelogfoldchangesTotal !=0,2] <- 1

x <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==0)]
y <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==0),1]

x1 <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==1)]
y1 <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==1),1]
plot(x,y,pch=16, xlab = "Index of genes", ylab = expression(paste("true"," ",log[2],"FC")),
     xlim = c(0,25000), ylim = c(-4,4),main = "Study 4")
points(x1,y1,col="red",pch = 16,ylim = c(-3,3),xlim = c(0,25000))








