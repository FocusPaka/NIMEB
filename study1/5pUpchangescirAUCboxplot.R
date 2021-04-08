
#单个样本DE ratio不同时的情况（循环）
rm(list=ls())
library('limma')
library(pROC)


source("DEmethods_noreplicate.R")
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

foldDiff <- 8
pUps <- c(0.5,0.7,0.9)
nc <- length(pUps)
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
        pUp <- pUps[j]
        xx <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                              foldDifference=foldDiff, pUp=pUp, 
                              pDifferential=0.6, empiricalDist=D[,1], 
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
        
        NIMEB_fdr[i,j] <- sum(pred_comm == FALSE)/sum(pred_all == FALSE)          #fdr
        NIMEB_pre[i,j] <- sum(pred_test==FALSE)/sum(pred_all==FALSE)              #precision
        NIMEB_sen[i,j] <- sum(pred_test==FALSE)/ndiff                             #sensitivity
        NIMEB_fs[i,j] <- 2*NIMEB_pre[i,j]*NIMEB_sen[i,j]/(NIMEB_pre[i,j]+NIMEB_sen[i,j])
        NIMEB_counts[i,j] <- sum(pred_all == FALSE)
        
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
        
        #---HTN method---#
        HTN_fdr[i,j] <- sum(res$pfull[comm,1] < 0.05)/sum(res$pfull[,1]<0.05)                              
        HTN_pre[i,j] <- sum(res$pfull[diff,1] < 0.05)/sum(res$pfull[,1] < 0.05)                            
        HTN_sen[i,j] <- sum(res$pfull[diff,1] < 0.05)/ndiff                                       
        HTN_fs[i,j] <- 2*HTN_pre[i,j]*HTN_sen[i,j]/(HTN_pre[i,j]+HTN_sen[i,j])
        HTN_counts[i,j] <- sum(res$pfull[,1]<0.05) 
        
        #-----HTN AUC---#
        check_ord_htn <- order(res$pfull[,1])
        roc_temp_htn <- roc(category[check_ord_htn], (-res$pfull[,1])[check_ord_htn])
        roc_obj_HTN[i,j] <- auc(roc_temp_htn)
        
        
        
        #---edgeR method---#
        edgeR_fdr[i,j] <- sum(res$pfull[comm,2] < 0.05)/sum(res$pfull[,2] < 0.05)
        edgeR_pre[i,j] <- sum(res$pfull[diff,2] < 0.05)/sum(res$pfull[,2] < 0.05)
        edgeR_sen[i,j] <- sum(res$pfull[diff,2] < 0.05)/ndiff
        edgeR_fs[i,j] <- 2*edgeR_pre[i,j]*edgeR_sen[i,j]/(edgeR_pre[i,j]+edgeR_sen[i,j])
        edgeR_counts[i,j] <- sum(res$pfull[,2] < 0.05)
        
        #----edgeR-like analysis AUC----#
        check_ord_edgeR <- order(res$pfull[,2])
        roc_temp_edgeR <- roc(category[check_ord_edgeR], (-res$pfull[,2])[check_ord_edgeR])
        roc_obj_edgeR[i,j] <- auc(roc_temp_edgeR)
        
        
        
        #---LibSize method---#
        LibSize_fdr[i,j] <- sum(res$pfull[comm,3] < 0.05)/sum(res$pfull[,3] < 0.05)
        LibSize_pre[i,j] <- sum(res$pfull[diff,3] < 0.05)/sum(res$pfull[,3] < 0.05)
        LibSize_sen[i,j] <- sum(res$pfull[diff,3] < 0.05)/ndiff
        LibSize_fs[i,j] <- 2*LibSize_pre[i,j]*LibSize_sen[i,j]/(LibSize_pre[i,j]+LibSize_sen[i,j])
        LibSize_counts[i,j] <- sum(res$pfull[,3] < 0.05)
        
        #------LibSize AUC------#
        check_ord_LibSize <- order(res$pfull[,3])
        roc_temp_LibSize <- roc(category[check_ord_LibSize], (-res$pfull[,3])[check_ord_LibSize])
        roc_obj_LibSize[i,j] <- auc(roc_temp_LibSize)
        
        
        #---DESeq method---#
        DESeq_fdr[i,j] <- sum(res$pfull[comm,4] < 0.05)/sum(res$pfull[,4] < 0.05)
        DESeq_pre[i,j] <- sum(res$pfull[diff,4] < 0.05)/sum(res$pfull[,4] < 0.05)
        DESeq_sen[i,j] <- sum(res$pfull[diff,4] < 0.05)/ndiff
        DESeq_fs[i,j] <- 2*DESeq_pre[i,j]*DESeq_sen[i,j]/(DESeq_pre[i,j]+DESeq_sen[i,j])
        DESeq_counts[i,j] <- sum(res$pfull[,4] < 0.05)
        
        #------DESeq AUC------#
        check_ord_DESeq <- order(res$pfull[,4])
        roc_temp_DESeq <- roc(category[check_ord_DESeq], (-res$pfull[,4])[check_ord_DESeq])
        roc_obj_DESeq[i,j] <- auc(roc_temp_DESeq)
        
        
        #---NOISeq method---#
        NOISeq_fdr[i,j] <- sum(res$pfull[comm,5] < 0.2)/sum(res$pfull[,5] < 0.2)
        NOISeq_pre[i,j] <- sum(res$pfull[diff,5] < 0.2)/sum(res$pfull[,5] < 0.2)
        NOISeq_sen[i,j] <- sum(res$pfull[diff,5] < 0.2)/ndiff
        NOISeq_fs[i,j] <- 2*NOISeq_pre[i,j]*NOISeq_sen[i,j]/(NOISeq_pre[i,j]+NOISeq_sen[i,j])
        NOISeq_counts[i,j] <- sum(res$pfull[,5] < 0.2)
        
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


###############################


##################################################################3
#AUC boxplot
auc_whole <- NULL
name <- NULL
for(i in 1:length(methods_name)){
    auc_whole <- rbind(auc_whole,get(paste("roc_obj_", methods_name[i],sep = "")))
    name <- c(name,rep(methods_name[i],reps))
}


data <- data.frame(matrix(NA,reps*length(methods_name),length(pUps)))
data[,1] <- name
data[,2] <- auc_whole[,1]
data[,3] <- auc_whole[,2]
data[,4] <- auc_whole[,3]
#data[,5] <- auc_whole[,4]

#data[,6] <- auc_whole[,5]
#data[,7] <- auc_whole[,6]

#data[,2:7] <- round(data[,2:7],2)

colnames(data) <- c("name","X1","X2","X3")   #,"X4")#,"X5","X6")
View(data)

require("RColorBrewer")
miscolores = brewer.pal(8,"Set3")[1:8]
par(mfrow=c(1,1))
oldpar=par(mar=c(5,3.6,1,0.2),mgp=c(1.7,0.5,0))


boxplot(X1~name,data,at=0:5, col = miscolores[1:6], xlim=c(0,19),ylim=c(0.5,1), ann = F,   cex.axis = 0.6,las=2)
boxplot(X2~name,data,at=7:12,add=T,cex.axis = 0.6,las=2,col = miscolores[1:6])
boxplot(X3~name,data,add=T,at=14:19,cex.axis = 0.6,las=2,col = miscolores[1:6])
#boxplot(X4~name,data,add=T,at=21:26,cex.axis = 0.6,las=2,col = miscolores[1:6])
title(xlab="Methods",ylab="AUC",line = 2.5)





abline(v=6,lty=3,col="pink")
abline(v=13,lty=3,col="pink")
#abline(v=20,lty=3,col="pink")

#abline(v=28,lty=3,col="pink")
#abline(v=35,lty=3,col="pink")




View(data)
library(ggplot2)
library(tidyr)
test <- data
test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
colnames(test) <- c("Methods","pUp=0.5","pUp=0.7","pUp=0.9")
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






