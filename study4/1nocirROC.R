#多个replicates=2,3,4,5,6,7的情况（循环）
rm(list=ls())
source("functions.R")
source("DEmethods_replicates.R")

#产生数据
data_primal <- read.table("Grimmond_lengths.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)

library('limma')

# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------



nreps <- c(2,2)
pUp <- 1
xx_new1 <- generateDataset2(commonTags=15000, uniqueTags=c(1000,800), 
                            empiricalDist = data_primal$EB,
                            lengthDist = data_primal$transcript_length,
                            pDifferential=0.6, foldDifference=4, pUp=pUp,
                            libLimits=c(.9,1.2)*1e6,nreps = nreps)

xx_new2 <- generateDataset2(commonTags = 10000, uniqueTags=c(2000,1000), 
                            libLimits=c(.9,1.2)*1e6,
                            empiricalDist = data_primal$EB,
                            lengthDist = data_primal$transcript_length, 
                            pDifferential = 0.4, pUp=.1, 
                            foldDifference= 8,nreps = nreps)


xx_new <- list()
xx_new$DATA <- rbind(xx_new1$DATA, xx_new2$DATA) 
xx_new$commonInd <- union(xx_new1$commonInd,xx_new2$commonInd+nrow(xx_new1$DATA))
xx_new$differentialInd <- union(xx_new1$differentialInd,xx_new2$differentialInd + 
                                    nrow(xx_new1$DATA))
xx_new$group <- xx_new1$group
xx_new$length <- c(xx_new1$length, xx_new2$length) 


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
condsAB <- xx$group
genelength <- c(xx_1$length,xx_2$length)

res <- runDEs(countsTable=countsTable, condsAB=condsAB, run_MEB=TRUE, train_id=id, 
              gamma=seq(1e-07, 2e-04, 5e-06), nu=0.01, reject_rate=0.1,
              run_Marioni0=FALSE, run_Marioni=FALSE,run_edgeR0=FALSE, 
              run_edgeR=TRUE,run_cloonan=FALSE,genelength=genelength, run_HTN=TRUE, 
              run_DESeq=TRUE, run_NOISeq=TRUE)


#----MEB method----#
pred_comm <- predict(res$MEBmodel,xx$DATA[comm,])           #无差异的基因
x_test <- xx$DATA[diff,]                                    #所有有差异的基因
pred_test <- predict(res$MEBmodel, x_test)
pred_all <- predict(res$MEBmodel, xx$DATA)                  #总的数据



#----MEB AUC-------#
check <- numeric(nrow(xx$DATA))
for(m in 1:nrow(xx$DATA)){
    check[m] <- decision_function(xx$DATA[m,], model=res$MEBmodel,gamma=res$gamma)-(res$MEBmodel)$rho
}
sum(check>0)           #no.TRUE   non-DE genes
summary(predict(res$MEBmodel,xx$DATA))

library(pROC)
check_ord_MEB <- order(check)
category <- numeric(nrow(xx$DATA))
category[diff] <- 1

roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

roc(category, check)



#################################################
#ROC curve
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))


plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE, asp=NA)
par(pty="s")




#---LibSize method---#
libsize.pvalues <- Poisson.model.new(countMatrix = countsTable, group1 = which(condsAB == 1), 
                                     group2 = which(condsAB == 2), calcFactor = FALSE)

libsize_roc <- roc(category, libsize.pvalues$stats$pval,smooth = T)
libsize_auc <- auc(libsize_roc)


plot(libsize_roc,col="green",add=T,legacy.axes=T)






#edgeR
roc_temp_edgeR <- roc(category, res$pfull[,4],smooth = T)
auc(roc_temp_edgeR)

plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)



#HTN
roc_temp_HTN <- roc(category, res$pfull[,6],smooth=T)
auc(roc_temp_HTN)

plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)




#DESeq
roc_temp_DESeq <- roc(category, res$pfull[,7],smooth=T)
auc(roc_temp_DESeq)

plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)



#NOISeq
roc_temp_NOISeq <- roc(category, res$pfull[,8],smooth=T)
auc(roc_temp_NOISeq)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T)


title("pUp=0.5")
legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESseq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)
















