#两个样本DE ratio不同时的情况  
rm(list=ls())

library('limma')
library(pROC)
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
xx_new1 <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                           foldDifference=4, pUp=0.9, 
                           pDifferential=0.6, empiricalDist=D[,1], 
                           libLimits=c(.9,1.2)*1e6)

xx_new2 <- generateDataset(commonTags=10000, uniqueTags=c(800,1500), 
                           foldDifference=8, pUp=0.1, 
                           pDifferential=0.3, empiricalDist=D[,1], 
                           libLimits=c(.9,1.2)*1e6)

xx_new1$trueFactors
xx_new2$trueFactors


xx_new <- list()
xx_new$DATA <- rbind(xx_new1$DATA, xx_new2$DATA) 
xx_new$commonInd <- union(xx_new1$commonInd,xx_new2$commonInd+nrow(xx_new1$DATA))
xx_new$differentialInd <- union(xx_new1$differentialInd,xx_new2$differentialInd + nrow(xx_new1$DATA))
xx_new$group <- xx_new1$group
xx_new$length <- rbind(xx_new1$length, xx_new2$length) 


k <- which(rowSums(xx_new$DATA) > 0 & rowMeans(xx_new$DATA) > 2)
xx <- takeSubset(xx_new, k)

k1 <- which(rowSums(xx_new1$DATA) > 0 & rowMeans(xx_new1$DATA) > 2)
xx_1 <- takeSubset(xx_new1, k1)

ci_1 <- xx_1$commonInd                                 #表达水平相同（包括成倍数的）
dii_1 <- intersect(ci_1,xx_1$differentialInd)          #表达水平成倍数的
comm_1 <- setdiff(ci_1,dii_1)                          #表达水平相同的
diff_1 <- xx_1$differentialInd                         #表达水平不同
ndiff_1 <- length(diff_1)
#################
k2 <- which(rowSums(xx_new2$DATA) > 0 & rowMeans(xx_new2$DATA) > 2)
xx_2 <- takeSubset(xx_new2, k2)

ci_2 <- xx_2$commonInd                                 #表达水平相同（包括成倍数的）
dii_2 <- intersect(ci_2,xx_2$differentialInd)          #表达水平成倍数的
comm_2 <- setdiff(ci_2,dii_2)                          #表达水平相同的
diff_2 <- xx_2$differentialInd                         #表达水平不同
ndiff_2 <- length(diff_2)
#########################################
ci <- xx$commonInd                               #表达水平相同（包括成倍数的）
dii <- intersect(ci,xx$differentialInd)          #表达水平成倍数的
comm <- setdiff(ci,dii)                          #表达水平相同的
diff <- xx$differentialInd                       #表达水平不同
ndiff <- length(diff)


id <- c(sample(comm_1, 500, replace=FALSE),sample(comm_2, 500, replace = FALSE)+nrow(xx_1$DATA))
x_train <- xx$DATA[id,]

countsTable <- xx$DATA

#calculate scaling factor
calcNormFactors_new(xx_new1$DATA,logratioTrim=0.3, sumTrim=0.05)
calcNormFactors_new(xx_new2$DATA,logratioTrim=0.3, sumTrim=0.05)
calcNormFactors_new(xx$DATA,logratioTrim=0.3, sumTrim=0.05)


colnames(xx$DATA) <- c('A','B')
condsAB <- colnames(xx$DATA)


res <- runDE(countsTable=countsTable, condsAB=condsAB, run_MEB=TRUE, 
             train_id=id, gamma=seq(0.0001,0.005,0.0001), nu=0.01, 
             reject_rate=0.1,run_HTN=TRUE, run_edgeR=TRUE, run_LibSize=TRUE, 
             run_DESeq=TRUE, run_NOISeq=TRUE)


#---MEB method---#
pred_comm <- predict(res$MEBmodel,xx$DATA[comm,])           #无差异的基因
summary(pred_comm)
x_test <- xx$DATA[diff,]                                    #所有有差异的基因
pred_test <- predict(res$MEBmodel, x_test)
summary(pred_test)
pred_all <- predict(res$MEBmodel, xx$DATA)                  #总的数据
summary(pred_all)



MEB_fdr <- sum(pred_comm == FALSE)/sum(pred_all == FALSE)          #fdr
MEB_fdr
MEB_pre <- sum(pred_test==FALSE)/sum(pred_all==FALSE)              #precision
MEB_sen <- sum(pred_test==FALSE)/ndiff                             #sensitivity
MEB_fs <- 2*MEB_pre*MEB_sen/(MEB_pre+MEB_sen)

###########################
pred_all_num <- numeric(length(pred_all))           #DE genes
pred_all_num[which(pred_all == TRUE)] <- 1          #non-DE genes
new_data_c <- cbind(xx$DATA, pred_all_num)



#----MEB AUC-------#
check <- numeric(nrow(xx$DATA))
for(m in 1:nrow(xx$DATA)){
    check[m] <- decision_function(xx$DATA[m,], model=res$MEBmodel,gamma=res$gamma)-(res$MEBmodel)$rho
}
#sum(check>0)           #no.TRUE   non-DE genes
check_ord_MEB <- order(check)
category <- numeric(nrow(xx$DATA))
category[diff] <- 1



#ROC curve
par(mfrow=c(1,1))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))

#MEB
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE,ylim = c(0,1))

par(pty="s")









#---HTN method---#
roc_temp_HTN <- roc(category, res$pfull[,1],smooth = T)
auc(roc_temp_HTN)

plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)





#---edgeR method---#
roc_temp_edgeR <- roc(category, res$pfull[,2],smooth = T)
auc(roc_temp_edgeR)

plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)




#---LibSize method---#
roc_temp_LibSize <- roc(category, res$pfull[,3],smooth = T)
auc(roc_temp_LibSize)

plot(roc_temp_LibSize,col="green",add=T,legacy.axes=T)






#---DESeq method---#
roc_temp_DESeq <- roc(category, res$pfull[,4],smooth = T)
auc(roc_temp_DESeq)

plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)




#---NOISeq method---#
roc_temp_NOISeq <- roc(category, res$pfull[,5],smooth = T)
auc(roc_temp_NOISeq)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T,print.auc=F)

#title(main = "d=4")

legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESeq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)








roc_obj_MEB
auc(roc_temp_HTN)
auc(roc_temp_edgeR)
auc(roc_temp_LibSize)
auc(roc_temp_DESeq)
auc(roc_temp_NOISeq)
