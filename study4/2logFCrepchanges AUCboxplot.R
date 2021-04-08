
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
replicates <- c(2,5,8)
nc <- length(replicates)
reps <- 10



fdr <- NULL
precision <- NULL
sensitivity <- NULL
f_score <- NULL
methods_name <- c("MEB","edgeR","LibrarySize", "HTN","DESeq","NOISeq")
for(i in 1:length(methods_name)){
    assign(paste(methods_name[i],"_fdr",sep = ""), matrix(0,reps,nc)) 
    assign(paste(methods_name[i], "_pre",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_sen",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_fs",sep = ""),matrix(0,reps,nc))
    assign(paste("roc_obj_", methods_name[i],sep = ""),matrix(0,reps,nc))
    
    fdr <- c(fdr,paste(methods_name[i],"_fdr",sep = ""))
    precision <- c(precision,paste(methods_name[i], "_pre",sep = ""))
    sensitivity <- c(sensitivity,paste(methods_name[i], "_sen",sep = ""))
    f_score <- c(f_score,paste(methods_name[i], "_fs",sep = ""))
    auc <- c(auc,paste("roc_obj_", methods_name[i],sep = ""))
}


for(i in 1:reps){
    for(j in 1:nc){
        nreps <- rep(replicates[j],2)
        xx_new1 <- generateDataset2(commonTags=15000, uniqueTags=c(1000,800), 
                                    empiricalDist = data_primal$EB,
                                    lengthDist = data_primal$transcript_length,
                                    pDifferential=0.6, foldDifference=4, pUp=0.9,
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
        xx_new$differentialInd <- union(xx_new1$differentialInd,xx_new2$differentialInd + nrow(xx_new1$DATA))
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
                      run_Marioni0=TRUE, run_Marioni=TRUE,run_edgeR0=TRUE, 
                      run_edgeR=TRUE,run_cloonan=TRUE,genelength=genelength, run_HTN=TRUE, 
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
        #sum(check>0)           #no.TRUE   non-DE genes
        check_ord_MEB <- order(check)
        category <- numeric(nrow(xx$DATA))
        category[diff] <- 1
        
        roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
        roc_obj_MEB[i,j] <- auc(roc_temp_MEB)
    
        #----edgeR-like analysis AUC----#
        check_ord_edgeR <- order(res$pfull[,4])
        roc_temp_edgeR <- roc(category[check_ord_edgeR], (-res$pfull[,4])[check_ord_edgeR])
        roc_obj_edgeR[i,j] <- auc(roc_temp_edgeR)
        
        
        #-----library size-------------#
        libsize.pvalues <- Poisson.model.new(countMatrix = countsTable, 
                                             group1 = which(condsAB == 1), 
                                             group2 = which(condsAB == 2), 
                                             calcFactor = FALSE)
        
        libsize_roc <- roc(category, libsize.pvalues$stats$pval)
        roc_obj_LibrarySize[i,j] <- auc(libsize_roc)

        #-----HTN AUC---#
        check_ord_htn <- order(res$pfull[,6])
        roc_temp_htn <- roc(category[check_ord_htn], (-res$pfull[,6])[check_ord_htn])
        roc_obj_HTN[i,j] <- auc(roc_temp_htn)
        
        
        #------DESeq AUC------#
        check_ord_DESeq <- order(res$pfull[,7])
        roc_temp_DESeq <- roc(category[check_ord_DESeq], (-res$pfull[,7])[check_ord_DESeq])
        roc_obj_DESeq[i,j] <- auc(roc_temp_DESeq)
        
        

        #------NOISeq  AUC------#
        check_ord_NOISeq <- order(-res$pfull[,8])
        roc_temp_NOISeq <- roc(category[check_ord_NOISeq], (res$pfull[,8])[check_ord_NOISeq])
        roc_obj_NOISeq[i,j] <- auc(roc_temp_NOISeq)
        
    }
}



###############################################################################
for(i in 1:length(methods_name)){
    print(fdr[i])
    print(get(paste(methods_name[i],"_fdr",sep = "")))
    
    print(precision[i])
    print(get(paste(methods_name[i],"_pre",sep = "")))
    
    print(sensitivity[i])
    print(get(paste(methods_name[i],"_sen",sep = "")))
    
    print(f_score[i])
    print(get(paste(methods_name[i],"_fs",sep = "")))
    
    print(auc[i])
    print(get(paste("roc_obj_", methods_name[i],sep = "")))
}


###############################
fdr_new <- NULL

pre_new <- NULL
sen_new <- NULL
fs_new <- NULL
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
    auc_new <- rbind(auc_new,apply(get(paste("roc_obj_", methods_name[i],sep = "")),2,mean))
}








################################################
#AUC boxplot
auc_whole <- NULL
name <- NULL
for(i in 1:length(methods_name)){
    auc_whole <- rbind(auc_whole,get(paste("roc_obj_", methods_name[i],sep = "")))
    name <- c(name,rep(methods_name[i],reps))
}

dim(auc_whole)
length(name)

name[1:10] <- rep("NIMEB", 10)

data <- data.frame(matrix(NA,reps*length(methods_name),length(replicates)))
data[,1] <- name
data[,2] <- auc_whole[,1]
data[,3] <- auc_whole[,2]
data[,4] <- auc_whole[,3]


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








View(data)
library(ggplot2)
library(tidyr)
test <- data
test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
colnames(test) <- c("Methods","Replicates number=2","Replicates number=5",
                    "Replicates number=8")
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

