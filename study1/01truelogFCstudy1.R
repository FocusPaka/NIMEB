#单个样本DE ratio不同时的情况  
rm(list=ls())

source("functions.R")


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
foldDiff <- 4
pUp <- 0.9
pDifferential <- 0.3

xx <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                      foldDifference=foldDiff, pUp=pUp, 
                      pDifferential=pDifferential, 
                      empiricalDist=D[,1], 
                      libLimits=c(.9,1.2)*1e6)
xx$trueFactors


#View(xx$truelog2foldchanges)
summary(xx$truelog2foldchanges)


truelogfoldchangesTotal <- xx$truelog2foldchanges

truelogfoldchangesTotalM <- matrix(0, nrow = length(truelogfoldchangesTotal),ncol = 2)

truelogfoldchangesTotalM[,1] <- truelogfoldchangesTotal
truelogfoldchangesTotalM[truelogfoldchangesTotal !=0,2] <- 1   #DE gene标志为1，non-DE标志为0

x <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==0)]    #先画non-DE
y <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==0),1]

x1 <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==1)]
y1 <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==1),1]
plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,15500), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 1")
points(x1,y1,col="red",pch = 16,ylim = c(-4,4),xlim = c(0,15500))



# library(ggplot2)
# View(truelogfoldchangesTotalM)
# test <- data.frame(index=1:nrow(truelogfoldchangesTotalM),
#                    y=truelogfoldchangesTotalM[,1],
#                    z=as.factor(truelogfoldchangesTotalM[,2]))
# 
# ggplot(data = test)+
#     geom_point(mapping = aes(x=test[,1],y=test[,2],
#                              color=test[,3]))

