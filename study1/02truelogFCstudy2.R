#两个样本DE ratio不同时的情况  
rm(list=ls())

library('limma')
#setwd("E:/paper/paper3_SVM/svm/code0510")
source("functions.R")
source("DEmethods_noreplicate.R")

# ------------------------------------
# table of counts from Marioni et al.
# ------------------------------------


load("LK_data.RData")

D <- as.matrix(MA.subsetA$M)
g <- as.character(MA.subsetA$genes$EnsemblGeneID)
o <- order(gsub("R[1-2]L[1-8]","",colnames(D)))


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------
pD <- 0.1


xx_new1 <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                           foldDifference=4, pUp=0.9, 
                           pDifferential=0.6, empiricalDist=D[,1], 
                           libLimits=c(.9,1.2)*1e6)

xx_new2 <- generateDataset(commonTags=10000, uniqueTags=c(800,1500), 
                           foldDifference=8, pUp=0.1, 
                           pDifferential=pD, empiricalDist=D[,1], 
                           libLimits=c(.9,1.2)*1e6)


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
plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,25000), xlab = "Index of genes", 
     ylab = expression(paste("true"," ",log[2],"FC")),main = "Study 2")
points(x1,y1,col="red",pch = 16,ylim = c(-4,4),xlim = c(0,25000))






