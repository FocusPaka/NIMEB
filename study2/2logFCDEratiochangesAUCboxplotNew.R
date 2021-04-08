#单个样本,2scaling factor, DE ratio不同时的情况（循环）
rm(list=ls())
library('limma')
source("functions.R")
source("DEmethods_noreplicate.R")

load("LK_data.RData")

library(pROC) 
# ------------------------------------
# table of counts from Marioni et al.
# ------------------------------------
D <- as.matrix(MA.subsetA$M)
g <- as.character(MA.subsetA$genes$EnsemblGeneID)
o <- order(gsub("R[1-2]L[1-8]","",colnames(D)))


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------
pDifferential <- c(0.3,0.5,0.7)
nc <- length(pDifferential)
reps <- 10

fdr <- NULL
precision <- NULL
sensitivity <- NULL
f_score <- NULL
auc_line <- NULL
counts_num <- NULL

methods_name <- c("MEB","HTN","edgeR","LibSize","DESeq","NOISeq")
for(i in 1:length(methods_name)){
    assign(paste(methods_name[i],"_fdr",sep = ""), matrix(0,reps,nc)) 
    assign(paste(methods_name[i], "_pre",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_sen",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_fs",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_AUC",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_counts",sep = ""),matrix(0,reps,nc))
    
    
    fdr <- c(fdr,paste(methods_name[i],"_fdr",sep = ""))
    precision <- c(precision,paste(methods_name[i], "_pre",sep = ""))
    sensitivity <- c(sensitivity,paste(methods_name[i], "_sen",sep = ""))
    f_score <- c(f_score,paste(methods_name[i], "_fs",sep = ""))
    auc_line <- c(auc_line,paste(methods_name[i], "_AUC",sep = ""))
    counts_num <- c(counts_num,paste(methods_name[i], "_counts",sep = ""))
}



for(i in 1:reps){
    for(j in 1:nc){
        pD <- pDifferential[j]
        xx_new1 <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                                   foldDifference=4, pUp=0.9, 
                                   pDifferential=pD, empiricalDist=D[,1], 
                                   libLimits=c(.9,1.2)*1e6)
        
        xx_new2 <- generateDataset(commonTags=10000, uniqueTags=c(800,1500), 
                                   foldDifference=8, pUp=0.1, 
                                   pDifferential=0.3, empiricalDist=D[,1], 
                                   libLimits=c(.9,1.2)*1e6)
        
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
        
        colnames(xx$DATA) <- c('A','B')
        condsAB <- colnames(xx$DATA)
        
        res <- runDE(countsTable=countsTable, condsAB=condsAB, run_MEB=TRUE, 
                     train_id=id, gamma=seq(0.0001,0.005,0.0001), nu=0.01, 
                     reject_rate=0.1,run_HTN=TRUE, run_edgeR=TRUE, run_LibSize=TRUE, 
                     run_DESeq=TRUE, run_NOISeq=TRUE)
        
        
        #---MEB method---#
        pred_comm <- predict(res$MEBmodel,xx$DATA[comm,])           #无差异的基因
        x_test <- xx$DATA[diff,]                                    #所有有差异的基因
        pred_test <- predict(res$MEBmodel, x_test)
        pred_all <- predict(res$MEBmodel, xx$DATA)                  #总的数据
        
        #---MEB AUC---#
        check <- numeric(nrow(xx$DATA))
        for(m in 1:nrow(xx$DATA)){
            check[m] <- decision_function(xx$DATA[m,], model=res$MEBmodel,gamma=res$gamma)-(res$MEBmodel)$rho
        }
        #sum(check>0)           #no.TRUE   non-DE genes
        check_ord_MEB <- order(check)
        category <- numeric(nrow(xx$DATA))
        category[diff] <- 1
        
        roc_temp_MEB <- roc(category, check)
        MEB_AUC[i,j] <- auc(roc_temp_MEB)
        MEB_counts[i,j] <- sum(pred_all == FALSE)
        
        
        #---HTN method---#
        #---HTN AUC---#
        roc_temp_HTN <- roc(category, res$pfull[,1])
        HTN_AUC[i,j] <- auc(roc_temp_HTN)
        HTN_counts[i,j] <- sum(res$pfull[,1]<0.05)
        
        #---edgeR AUC---#
        roc_temp_edgeR <- roc(category, res$pfull[,2])
        edgeR_AUC[i,j] <- auc(roc_temp_edgeR)
        edgeR_counts[i,j] <- sum(res$pfull[,2] < 0.05)
        
        #---LibSize AUC---#
        roc_temp_LibSize <- roc(category, res$pfull[,3])
        LibSize_AUC[i,j] <- auc(roc_temp_LibSize)
        LibSize_counts[i,j] <- sum(res$pfull[,3] < 0.05)
        
        
        #---DESeq AUC---#
        roc_temp_DESeq <- roc(category, res$pfull[,4])
        DESeq_AUC[i,j] <- auc(roc_temp_DESeq)
        DESeq_counts[i,j] <- sum(res$pfull[,4] < 0.05)
        
        #---NOISeq AUC---#
        roc_temp_NOISeq <- roc(category, res$pfull[,5])
        NOISeq_AUC[i,j] <- auc(roc_temp_NOISeq)
        NOISeq_counts[i,j] <- sum(res$pfull[,5] < 0.2)
        
    }
}




############################################################
for(i in 1:length(methods_name)){
    print(fdr[i])
    print(get(paste(methods_name[i],"_fdr",sep = "")))
    
    print(precision[i])
    print(get(paste(methods_name[i],"_pre",sep = "")))
    
    print(sensitivity[i])
    print(get(paste(methods_name[i],"_sen",sep = "")))
    
    print(f_score[i])
    print(get(paste(methods_name[i],"_fs",sep = "")))
    
    print(auc_line[i])
    print(get(paste(methods_name[i],"_AUC",sep = "")))
    
    print(counts_num[i])
    print(get(paste(methods_name[i],"_counts",sep = "")))
}


###############################
fdr_new <- NULL
pre_new <- NULL
sen_new <- NULL
fs_new <- NULL
auc_new <- NULL
counts_num_new <- NULL
for(i in 1:length(methods_name)){
    fdr_new <- rbind(fdr_new,apply(get(paste(methods_name[i],"_fdr",sep = "")),
                                   2,mean))
    pre_new <- rbind(pre_new,apply(get(paste(methods_name[i],"_pre",sep = "")),
                                   2,mean))
    sen_new <- rbind(sen_new,apply(get(paste(methods_name[i],"_sen",sep = "")),
                                   2,mean))
    fs_new <- rbind(fs_new,apply(get(paste(methods_name[i],"_fs",sep = "")),
                                 2,mean))
    
    auc_new <- rbind(auc_new,apply(get(paste(methods_name[i],"_AUC",sep = "")),
                                   2,mean))
    counts_num_new <- rbind(counts_num_new,apply(get(paste(methods_name[i],"_counts",sep = "")),
                                                 2,mean))
}



##############################################
#AUC boxplot   只要一部分
auc_whole <- NULL
name <- NULL
for(i in 1:length(methods_name)){
    auc_whole <- rbind(auc_whole,get(paste(methods_name[i],"_AUC",sep = "")))
    name <- c(name,rep(methods_name[i],reps))
}

name[1:10] <- rep("NIMEB", 10)
#name[61:80] <- rep("Library Size", 20)


data <- data.frame(matrix(NA,reps*length(methods_name),length(pDifferential)))
data[,1] <- name
data[,2] <- auc_whole[,1]
data[,3] <- auc_whole[,2]
data[,4] <- auc_whole[,3]

#data[,2:4] <- round(data[,2:4],2)

colnames(data) <- c("name","X1","X2","X3")


require("RColorBrewer")
miscolores = brewer.pal(8,"Set3")[1:8]
par(mfrow=c(1,1))
oldpar=par(mar=c(5,4,1,0.6),mgp=c(1.7,0.5,0))
boxplot(X1~name,data,at=0:5,xlim=c(0,19),ylim = c(0,1),ann = F, col = miscolores[1:6], cex.axis = 0.6,las=2)
boxplot(X2~name,data,at=7:12,add=T,cex.axis = 0.6,las=2,col = miscolores[1:6])
boxplot(X3~name,data,add=T,at=14:19,cex.axis = 0.6,las=2,col = miscolores[1:6])
title(xlab="Methods",ylab="AUC",line = 3.1)

abline(v=6,lty=3,col="pink")
abline(v=13,lty=3,col="pink")
#abline(v=21,lty=3,col="pink")
#abline(v=28,lty=3,col="pink")
#abline(v=35,lty=3,col="pink")























