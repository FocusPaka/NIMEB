
#单个样本DE ratio不同时的情况（循环）
rm(list=ls())
library('limma')
library(pROC)
#setwd("E:/paper/paper3_SVM/svm/code0510")
#source("data_analysis.R")
source("DEmethods_noreplicate.R")
source("functions.R")
# ------------------------------------
# table of counts from Marioni et al.
# ------------------------------------
#path <- "E:/paper/paper3_SVM/svm/code0510/data/"
load("LK_data.RData")
D <- as.matrix(MA.subsetA$M)
g <- as.character(MA.subsetA$genes$EnsemblGeneID)
o <- order(gsub("R[1-2]L[1-8]","",colnames(D)))


# ------------------------------------
# simulate data from empirical distribution of counts
# ------------------------------------

foldDiff <- 4
pUp <- 0.9
pDifferential <- c(0.3,0.5,0.7)
nc <- length(pDifferential)
reps <- 10

fdr <- NULL
precision <- NULL
sensitivity <- NULL
f_score <- NULL
counts_num <- NULL
auc <- NULL
methods_name <- c("NIMEB","HTN","edgeR","LibSize","DESeq","NOISeq")
for(i in 1:length(methods_name)){
    assign(paste(methods_name[i], "_fdr",sep = ""), matrix(0,reps,nc)) 
    assign(paste(methods_name[i], "_pre",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_sen",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_fs",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_counts",sep = ""),matrix(0,reps,nc))
    assign(paste("roc_obj_", methods_name[i],sep = ""),matrix(0,reps,nc))
    
    
    fdr <- c(fdr,paste(methods_name[i],"_fdr",sep = ""))
    precision <- c(precision,paste(methods_name[i], "_pre",sep = ""))
    sensitivity <- c(sensitivity,paste(methods_name[i], "_sen",sep = ""))
    f_score <- c(f_score,paste(methods_name[i], "_fs",sep = ""))
    counts_num <- c(counts_num,paste(methods_name[i], "_counts",sep = ""))
    auc <- c(auc,paste("roc_obj_", methods_name[i],sep = ""))
    
}



for(i in 1:reps){
    for(j in 1:nc){
        pD <- pDifferential[j]
        xx <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                              foldDifference=foldDiff, pUp=pUp, 
                              pDifferential=pD, empiricalDist=D[,1], 
                              libLimits=c(.9,1.2)*1e6)
        
        k <- which(rowSums(xx$DATA) > 0 & rowMeans(xx$DATA) > 2)
        xx <- takeSubset(xx, k)
        
        ci <- xx$commonInd                               #表达水平相同（包括成倍数的）
        dii <- intersect(ci,xx$differentialInd)          #表达水平成倍数的
        comm <- setdiff(ci, dii)                         #表达水平相同的
        diff <- xx$differentialInd                       #表达水平不同
        ndiff <- length(diff)
        
        id <- sample(comm, 1000, replace=FALSE)
        countsTable <- xx$DATA
        condsAB <- c('A','B')
        
        res <- runDE(countsTable=countsTable, condsAB=condsAB, run_MEB=TRUE, 
                     train_id=id, gamma=seq(0.0001,0.005,0.0001), nu=0.01, 
                     reject_rate=0.05,run_HTN=TRUE, run_edgeR=TRUE, run_LibSize=TRUE, 
                     run_DESeq=TRUE, run_NOISeq=TRUE)
        
        
        #---MEB method---#
        pred_comm <- predict(res$MEBmodel,xx$DATA[comm,])           #无差异的基因
        x_test <- xx$DATA[diff,]                                    #所有有差异的基因
        pred_test <- predict(res$MEBmodel, x_test)
        pred_all <- predict(res$MEBmodel, xx$DATA)                  #总的数据
        
        
        #----MEB AUC-------#
        check <- numeric(nrow(xx$DATA))
        for(m in 1:nrow(xx$DATA)){
            check[m] <- decision_function(xx$DATA[m,], model=res$MEBmodel,gamma=res$gamma)-(res$MEBmodel)$rho
        }
        #sum(check>0)           #no.TRUE   non-DE genes
        check_ord_MEB <- order(check)
        category <- c(rep(1,length(dii)),rep(0,length(comm)),rep(1,nrow(xx$DATA)-length(ci)))
        roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
        roc_obj_NIMEB[i,j] <- auc(roc_temp_MEB)
        
        
        
        #-----HTN AUC---#
        check_ord_htn <- order(res$pfull[,1])
        roc_temp_htn <- roc(category[check_ord_htn], (-res$pfull[,1])[check_ord_htn])
        roc_obj_HTN[i,j] <- auc(roc_temp_htn)
        
        
        #----edgeR-like analysis AUC----#
        check_ord_edgeR <- order(res$pfull[,2])
        roc_temp_edgeR <- roc(category[check_ord_edgeR], (-res$pfull[,2])[check_ord_edgeR])
        roc_obj_edgeR[i,j] <- auc(roc_temp_edgeR)
        
        
        
        #------LibSize AUC------#
        check_ord_LibSize <- order(res$pfull[,3])
        roc_temp_LibSize <- roc(category[check_ord_LibSize], (-res$pfull[,3])[check_ord_LibSize])
        roc_obj_LibSize[i,j] <- auc(roc_temp_LibSize)
        
        
        #------DESeq AUC------#
        check_ord_DESeq <- order(res$pfull[,4])
        roc_temp_DESeq <- roc(category[check_ord_DESeq], (-res$pfull[,4])[check_ord_DESeq])
        roc_obj_DESeq[i,j] <- auc(roc_temp_DESeq)
        
        
        #------NOISeq  AUC------#
        check_ord_NOISeq <- order(-res$pfull[,5])
        roc_temp_NOISeq <- roc(category[check_ord_NOISeq], (res$pfull[,5])[check_ord_NOISeq])
        roc_obj_NOISeq[i,j] <- auc(roc_temp_NOISeq)
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
    
    print(counts_num[i])
    print(get(paste(methods_name[i],"_counts",sep = "")))
    
    print(auc[i])
    print(get(paste("roc_obj_", methods_name[i],sep = "")))
}


###############################
fdr_new <- NULL
pre_new <- NULL
sen_new <- NULL
fs_new <- NULL
counts_num_new <- NULL
auc_new <- NULL

for(i in 1:length(methods_name)){
    fdr_new <- rbind(fdr_new,apply(get(paste(methods_name[i],"_fdr",sep = "")),
                                   2,mean))
    pre_new <- rbind(pre_new,apply(get(paste(methods_name[i],"_pre",sep = "")),
                                   2,mean))
    sen_new <- rbind(sen_new,apply(get(paste(methods_name[i],"_sen",sep = "")),
                                   2,mean))
    fs_new <- rbind(fs_new,apply(get(paste(methods_name[i],"_fs",sep = "")),
                                 2,mean))
    counts_num_new <- rbind(counts_num_new,apply(get(paste(methods_name[i],"_counts",sep = "")),
                                                 2,mean))
    auc_new <- rbind(auc_new,apply(get(paste("roc_obj_", methods_name[i],sep = "")),2,mean))
    
}



##################################################################3
#AUC boxplot
auc_whole <- NULL
name <- NULL
for(i in 1:length(methods_name)){
    auc_whole <- rbind(auc_whole,get(paste("roc_obj_", methods_name[i],sep = "")))
    name <- c(name,rep(methods_name[i],reps))
}


#name[61:80] <- rep("Library Size", 20)
#name[1:20] <- rep("NIMEB", 20)

data <- data.frame(matrix(NA,reps*length(methods_name),length(pDifferential)))
data[,1] <- name
data[,2] <- auc_whole[,1]
data[,3] <- auc_whole[,2]
data[,4] <- auc_whole[,3]
#data[,5] <- auc_whole[,4]

#data[,6] <- auc_whole[,5]
#data[,7] <- auc_whole[,6]

#data[,2:7] <- round(data[,2:7],2)

colnames(data) <- c("name","X1","X2","X3")

require("RColorBrewer")
miscolores = brewer.pal(8,"Set3")[1:8]
par(mfrow=c(1,1))
oldpar=par(mar=c(5,3.6,1,0.2),mgp=c(1.7,0.5,0))


boxplot(X1~name,data,at=0:5, col = miscolores[1:6], xlim=c(0,19),ylim=c(0,1), ann = F,   cex.axis = 0.6,las=2)
boxplot(X2~name,data,at=7:12,add=T,cex.axis = 0.6,las=2,col = miscolores[1:6])
boxplot(X3~name,data,add=T,at=14:19,cex.axis = 0.6,las=2,col = miscolores[1:6])
#boxplot(X4~name,data,add=T,at=21:26,cex.axis = 0.6,las=2,col = miscolores[1:6])
title(xlab="Methods",ylab="AUC",line = 2.5)



abline(v=6,lty=3,lwd=2,col="pink")
abline(v=13,lty=3,lwd=2,col="pink")









library(ggplot2)
library(tidyr)
View(data)
test <- data
test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
View(test)


colnames(test) <- c("Methods","DE rate=0.3","DE rate=0.5","DE rate=0.7")
#colnames(test) <- c("Methods","pUp=0.5","pUp=0.7","pUp=0.9")

#View(test)
#扁变长
test_gather <- gather(data=test,
                      key=pUp,
                      value=AUC,
                      -Methods)
#View(test_gather)
ggplot(data=test_gather)+
    geom_boxplot(mapping=aes(x=Methods,y=AUC,color=Methods))+
    facet_wrap(~pUp)+
    coord_cartesian(ylim = c(0, 1))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

