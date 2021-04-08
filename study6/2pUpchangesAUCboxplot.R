rm(list = ls())


library('pROC')
library('compcodeR')
library('edgeR')
library('DESeq2')
library('e1071')
library('NOISeq')
source("functions.R")
source("DEmethods.R")


#生成第一个数据集
n.vars1 <- 15000
samples.pc.comm <- 2
samples.pc1 <- samples.pc.comm
n.diffexp1 <- n.vars1*0.6
fraction.upregulated1 <- 1
effect.size1 <- 4


n.vars2 <- 10000
samples.pc2 <- samples.pc.comm
n.diffexp2 <- n.vars2*0.4
fraction.upregulated.vector <- c(0.1,0.3,0.5)
fraction.upregulated2 <- fraction.upregulated.vector
effect.size2 <- 8



reps <- 10
nc <- length(fraction.upregulated.vector)



fdr <- NULL
precision <- NULL
sensitivity <- NULL
f_score <- NULL
auc <- NULL
methods_name <- c("NIMEB","libsize", "edgeR", "HTN","DESeq","NOISeq")

for (i in 1:length(methods_name)) {
    assign(paste(methods_name[i],"_fdr",sep = ""), matrix(0,reps,nc)) 
    assign(paste(methods_name[i], "_pre",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_sen",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_fs",sep = ""),matrix(0,reps,nc))
    assign(paste(methods_name[i], "_auc",sep = ""),matrix(0,reps,nc))
    
    fdr <- c(fdr,paste(methods_name[i],"_fdr",sep = ""))
    precision <- c(precision,paste(methods_name[i], "_pre",sep = ""))
    sensitivity <- c(sensitivity,paste(methods_name[i], "_sen",sep = ""))
    f_score <- c(f_score,paste(methods_name[i], "_fs",sep = ""))
    auc <- c(auc,paste(methods_name[i], "_auc", sep = ""))
}






for (i in 1:reps) {
    for (j in 1:nc) {
        B1 <- generateSyntheticData(dataset = "B1", n.vars = n.vars1, 
                                    samples.per.cond = samples.pc1, 
                                    n.diffexp = n.diffexp1, 
                                    repl.id = 1, seqdepth = 1e6,
                                    minfact = 0.7, maxfact = 1.4,
                                    relmeans = "auto",
                                    dispersions = "auto",
                                    fraction.upregulated = fraction.upregulated1, 
                                    between.group.diffdisp = FALSE, 
                                    filter.threshold.total = 1, 
                                    filter.threshold.mediancpm = 0, 
                                    fraction.non.overdispersed = 1, 
                                    random.outlier.high.prob = 0,
                                    random.outlier.low.prob = 0,
                                    single.outlier.high.prob = 0,
                                    single.outlier.low.prob = 0,
                                    effect.size = effect.size1,
                                    output.file = NULL)
        
        
        
        ndiff1 <- sum(B1@variable.annotations$differential.expression == 1)    #第一个数据集中差异表达基因的个数
        ndata1 <- nrow(B1@count.matrix)      #第一个数据集的总个数
        
        
        commDE1 <- B1@count.matrix[ndiff1:ndata1,]    #第一个数据集中非差异表达基因的数据
        diffDE1 <- B1@count.matrix[1:ndiff1,]         #第一个数据集中差异表达基因的数据
        dat1 <- B1@count.matrix                        #第一个数据集数据
        
        
        
        #------------------------------------------------------------------------------
        B2 <- generateSyntheticData(dataset = "B2", n.vars = n.vars2, 
                                    samples.per.cond = samples.pc2, 
                                    n.diffexp = n.diffexp2, 
                                    repl.id = 1, seqdepth = 1e6,
                                    minfact = 0.7, maxfact = 1.4,
                                    relmeans = "auto",
                                    dispersions = "auto",
                                    fraction.upregulated = fraction.upregulated2[j], 
                                    between.group.diffdisp = FALSE, 
                                    filter.threshold.total = 1, 
                                    filter.threshold.mediancpm = 0, 
                                    fraction.non.overdispersed = 1, 
                                    random.outlier.high.prob = 0,
                                    random.outlier.low.prob = 0,
                                    single.outlier.high.prob = 0,
                                    single.outlier.low.prob = 0,
                                    effect.size = effect.size2,
                                    output.file = NULL)
        
        
        ndiff2 <- sum(B2@variable.annotations$differential.expression == 1)      #第二个数据集中差异表达基因的个数
        ndata2 <- nrow(B2@count.matrix)                #第二个数据集的总数
        
        
        commDE2 <- B2@count.matrix[ndiff2:ndata2,]
        diffDE2 <- B2@count.matrix[1:ndiff2,]
        dat2 <- B2@count.matrix
        
        
        
        #-----------------------------NIMEB--------------------------------------------
        id1 <- sample((ndiff1+1):ndata1, 500, replace = FALSE)
        id2 <- sample((ndiff2+1):ndata2, 500, replace = FALSE)
        x_train <- rbind(dat1[id1,],dat2[id2,])
        
        dat <- rbind(dat1, dat2)
        condsAB <- B1@sample.annotations$condition
        id <- c(id1,id2+nrow(dat1))
        
        res <- runDEs(countsTable = dat, condsAB = condsAB, run_NIMEB = TRUE, 
                      train_id = id, gamma = seq(0,1e-2,1e-5),
                      run_librarySize = TRUE, run_edgeR = TRUE, run_HTN = TRUE, 
                      run_DESeq2 = TRUE, run_NOISeq = TRUE)
        
        
        
        #-----------------------------NIMEB--------------------------------------------
        category <- c(B1@variable.annotations$differential.expression, 
                      B2@variable.annotations$differential.expression)
        
        #----NIMEB AUC-------#
        check <- numeric(nrow(dat))
        for (k in 1:nrow(dat)) {
            check[k] <- decision_function(dat[k,], model=res$model,gamma=res$model$gamma) - res$model$rho
        }
        
        check_ord_MEB <- order(check)
        NIMEB_roc <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
        NIMEB_auc[i,j] <- auc(NIMEB_roc)
        
        
        
        #---------------------------Library size---------------------------------------
        libsize_roc <- roc(category, res$pfull[,1],smooth = T)
        libsize_auc[i,j] <- auc(libsize_roc)
        
        
        
        #------------------------------edgeR--------------------------------------------
        edgeR_roc <- roc(category, res$pfull[,2],smooth = T)
        edgeR_auc[i,j] <- auc(edgeR_roc)
        
        
        
        #----------------------------HTN--------------------------------------------------#
        HTN_roc <- roc(category, res$pfull[,3],smooth = T)
        HTN_auc[i,j] <- auc(HTN_roc)
        
        
        
        #------------------------------DESeq2------------------------------------------
        
        DESeq_roc <- roc(category, res$pfull[,4],smooth = T)
        DESeq_auc[i,j] <- auc(DESeq_roc)
        
        
        #------------------------------NOIseq------------------------------------------
        roc_NOISeq <- roc(category, res$pfull[,5],smooth = T)
        NOISeq_auc[i,j] <- auc(roc_NOISeq)
        
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
    print(get(paste(methods_name[i], "_auc",sep = "")))
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
    auc_new <- rbind(auc_new,apply(get(paste(methods_name[i], "_auc",sep = "")),2,mean))
}



#######################################
#可视化
#AUC boxplot
auc_whole <- NULL
name <- NULL
for(i in 1:length(methods_name)){
    auc_whole <- rbind(auc_whole,get(paste(methods_name[i],"_auc",sep = "")))
    name <- c(name,rep(methods_name[i],reps))
}




data <- data.frame(matrix(NA,reps*length(methods_name),nc))
data[,1] <- name
data[,2] <- auc_whole[,1]
data[,3] <- auc_whole[,2]
data[,4] <- auc_whole[,3]
data[,2:4] <- round(data[,2:4],2)

colnames(data) <- c("name","X1","X2","X3")



oldpar=par(mar=c(4,3.6,1,0.1),mgp=c(1.7,0.5,0))
require("RColorBrewer")
miscolores = brewer.pal(8,"Set3")[1:8]
par(mfrow=c(1,1))
#oldpar=par(mar=c(5,3.6,1,0.2),mgp=c(1.7,0.5,0))



boxplot(X1~name,data,at=0:5,xlim=c(0,19),ylim=c(0,1),
        ann = F,col = c(miscolores[1:6]), cex.axis = 0.6,las=2)

boxplot(X2~name,data,at=7:12,add=T,cex.axis = 0.6,las=2,
        col = c(miscolores[1:6]))
boxplot(X3~name,data,add=T,at=14:19,cex.axis = 0.6,las=2,
        col = c(miscolores[1:6]))
title(xlab="Methods",ylab="AUC",line = 2.5)

abline(v=6,lty=3,lwd=2,col="pink")
abline(v=13,lty=3,lwd=2,col="pink")






View(data)
library(ggplot2)
library(tidyr)
test <- data
test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
colnames(test) <- c("Methods","pUp=0.1","pUp=0.3","pUp=0.5")

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





