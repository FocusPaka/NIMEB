rm(list = ls())


library('pROC')
library('compcodeR')
library('edgeR')
library('DESeq2')
library('e1071')
library('NOISeq')
source("functions.R")
source("DEmethods.R")


n.vars <- 15000
samples.pc <- c(2,5,8)
pD <- 0.6
n.diffexp <- n.vars * pD
fraction.upregulated <- 1
effect.size <- 8
reps <- 10
dispersions <- rep(0.01,n.vars)
nc <- length(samples.pc)



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
        B_625_625 <- generateSyntheticData(dataset = paste0("B", i, j), 
                                           n.vars = n.vars, 
                                           samples.per.cond = samples.pc[j], 
                                           n.diffexp = n.diffexp, 
                                           repl.id = 1, seqdepth = 1e6,
                                           minfact = 0.7, maxfact = 1.4,
                                           relmeans = "auto",
                                           dispersions = dispersions,
                                           fraction.upregulated = fraction.upregulated, 
                                           between.group.diffdisp = FALSE, 
                                           filter.threshold.total = 1, 
                                           filter.threshold.mediancpm = 0, 
                                           fraction.non.overdispersed = 0, 
                                           random.outlier.high.prob = 0,
                                           random.outlier.low.prob = 0,
                                           single.outlier.high.prob = 0,
                                           single.outlier.low.prob = 0,
                                           effect.size = effect.size,
                                           output.file = NULL)
        
        
        ndiff <- sum(B_625_625@variable.annotations$differential.expression == 1)
        ndata <- nrow(B_625_625@count.matrix)
        category <- B_625_625@variable.annotations$differential.expression
        
        
        commDE <- B_625_625@count.matrix[(ndiff + 1):ndata,]      #non-DE genes
        diffDE <- B_625_625@count.matrix[1:ndiff,]          #DE genes
        dat <- B_625_625@count.matrix                       #data
        
        condsAB <- B_625_625@sample.annotations$condition
        print(head(dat))
        
        id <- sample((ndiff + 1):ndata, 1000, replace = FALSE) 
        
        res <- runDEs(countsTable = dat, condsAB = condsAB, run_NIMEB = TRUE, 
                      train_id = id, gamma = seq(0,1e-3,1e-6),
                      run_librarySize = TRUE, run_edgeR = TRUE, run_HTN = TRUE, 
                      run_DESeq2 = TRUE, run_NOISeq = TRUE)
        
        
        
        #-----------------------------NIMEB--------------------------------------------
        check <- numeric(nrow(dat))
        for (k in 1:nrow(dat)) {
            check[k] <- decision_function(dat[k,], model=res$model,gamma=res$model$gamma) - res$model$rho
        }
        
        check_ord_MEB <- order(check)
        NIMEB_roc <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth=T)
        NIMEB_auc[i,j] <- auc(NIMEB_roc)
        
        
        
        #---------------------------Library size---------------------------------------
        libsize_roc <- roc(category, res$pfull[,1],smooth=T)
        libsize_auc[i,j] <- auc(libsize_roc)
        
        
        
        #------------------------------edgeR--------------------------------------------

        edgeR_roc <- roc(category, res$pfull[,2],smooth=T)
        edgeR_auc[i,j] <- auc(edgeR_roc)
        
        
        
        #----------------------------HTN--------------------------------------------------#
        
        HTN_roc <- roc(category, res$pfull[,3],smooth=T)
        HTN_auc[i,j] <- auc(HTN_roc)
        
        
        
        #------------------------------DESeq2------------------------------------------
        
        
        DESeq_roc <- roc(category, res$pfull[,4],smooth=T)
        DESeq_auc[i,j] <- auc(DESeq_roc)
        
        
        #------------------------------NOIseq------------------------------------------
          
        roc_NOISeq <- roc(category, res$pfull[,5],smooth=T)
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
    auc_new <- rbind(auc_new,apply(get(paste(methods_name[i], "_auc",sep = "")),
                                   2,mean))
}






##################################################################3
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



