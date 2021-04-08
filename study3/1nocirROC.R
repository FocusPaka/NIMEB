#多个replicates的情况
rm(list=ls())


library(pROC)

source("functions.R")
source("DEmethods_replicates.R")

#产生数据
data_primal <- read.table("Grimmond_lengths.txt",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------
pD <- 0.6
nreps <- c(2,2)
xx <- generateDataset2(commonTags = 15000, uniqueTags=c(1000,800),
                       group=c(1,2), libLimits=c(.9,1.2)*1e6,
                       empiricalDist = data_primal$EB,
                       lengthDist = data_primal$transcript_length,
                       pDifferential = pD, pUp=0.9,
                       foldDifference= 4, nreps= nreps)


k <- which(rowSums(xx$DATA) > 0 & rowMeans(xx$DATA) > 2)
xx <- takeSubset(xx, k)

ci <- xx$commonInd                               #表达水平相同（包括成倍数的）
dii <- intersect(ci,xx$differentialInd)          #表达水平成倍数的
comm <- setdiff(ci,dii)                          #表达水平相同的
diff <- xx$differentialInd                       #表达水平不同
ndiff <- length(diff)

id <- sample(comm, 1000, replace=FALSE)
x_train <- xx$DATA[id,]
countsTable <- xx$DATA
condsAB <- xx$group
genelength <- xx$length


res <- runDEs(countsTable=countsTable, condsAB=condsAB, run_MEB=TRUE, train_id=id, 
              gamma=seq(1e-07, 2e-04, 5e-06), nu=0.01, reject_rate=0.05,
              run_Marioni0=FALSE, 
              run_Marioni=FALSE,
              run_edgeR0=FALSE, 
              run_edgeR=TRUE,run_cloonan=FALSE,genelength=genelength, 
              run_HTN=TRUE, 
              run_DESeq=TRUE, run_NOISeq=TRUE)



#################################################################
#----MEB AUC-------#
check <- numeric(nrow(xx$DATA))
for(m in 1:nrow(xx$DATA)){
    check[m] <- decision_function(xx$DATA[m,], model=res$MEBmodel,gamma=res$gamma)-(res$MEBmodel)$rho
}
sum(check>0)           #no.TRUE   non-DE genes
summary(predict(res$MEBmodel,xx$DATA))


check_ord_MEB <- order(check)
category <- c(rep(1,length(dii)),rep(0,length(comm)),rep(1,nrow(xx$DATA)-length(ci)))   #DE genes equal to 1
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

roc(category, check)


###########################################################
par(mfrow=c(1,1))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))

#MEB
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

#plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE, asp=NA, ann = F)
plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE,lwd=3, xlim=c(1,0),ylim = c(0,1))

par(pty="s")


#HTN
roc_temp_HTN <- roc(category, res$pfull[,6],smooth = T)
auc(roc_temp_HTN)

plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)



#---------------------------Library size---------------------------------------
libsize.pvalues <- Poisson.model.new(countMatrix = countsTable, group1 = which(condsAB == 1), 
                                     group2 = which(condsAB == 2), calcFactor = FALSE)

libsize_roc <- roc(category, libsize.pvalues$stats$pval)
libsize_auc <- auc(libsize_roc)


plot(libsize_roc,col="green",add=T,legacy.axes=T)





#edgeR
roc_temp_edgeR <- roc(category, res$pfull[,4],smooth = T)
auc(roc_temp_edgeR)

plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)






#DESeq
roc_temp_DESeq <- roc(category, res$pfull[,7],smooth = T)
auc(roc_temp_DESeq)

plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)



#NOISeq
roc_temp_NOISeq <- roc(category, res$pfull[,8],smooth = T)
auc(roc_temp_NOISeq)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T)


legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESseq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)








