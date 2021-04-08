#单个样本DE ratio不同时的情况  
rm(list=ls())
library('e1071')
library('limma')
library('edgeR')
library('statmod')
library('DESeq')
library('NOISeq')
library('pROC')



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
pUp <- 0.9
pDifferential <- 0.55

xx <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                      foldDifference=foldDiff, pUp=pUp, 
                      pDifferential=pDifferential, empiricalDist=D[,1], 
                      libLimits=c(.9,1.2)*1e6)

k <- which(rowSums(xx$DATA) > 0 & rowMeans(xx$DATA) > 2)
xx <- takeSubset(xx, k)


ci <- xx$commonInd                               #表达水平相同（包括成倍数的）
dii <- intersect(ci,xx$differentialInd)          #表达水平成倍数的
comm <- setdiff(ci, dii)                         #表达水平相同的
diff <- xx$differentialInd                       #表达水平不同
ndiff <- length(diff)


############################################
id <- sample(comm, 1000, replace=FALSE)
x_train <- xx$DATA[id,]
#gamma <- seq(0,0.01,0.0001)

#gamma <- seq(0,0.01,length.out = 10000)

gamma <- seq(0,0.01,length.out = 10000)





test_error <- numeric(length(gamma))


train_error <- numeric(length(gamma))
for(i in 1:length(gamma)){
    model <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification",
                 kernel = "radial", gamma = gamma[i],
                 nu = 0.01, tolerance = 0.001,
                 shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                 na.action = na.omit)
    
    train_error[i] <- 1 - model$tot.accuracy/100
}


#######################################################
#gamma_num_new <- which.min(abs(train_error-0.05))
#train_error[gamma_num_new]

for(i in 1:length(gamma)){   
    model_new <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification", 
                     kernel = "radial", gamma = gamma[i],
                     nu = 0.01, tolerance = 0.001, 
                     shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                     na.action = na.omit)
    
    # summary(model)
    # pred_train <- predict(model, x_train)
    # summary(pred_train) 
    # train_error_all <- 1-sum(pred_train)/nrow(x_train)
    # train_error_all
    # 
    # ###########################################
    # #评估
    pred_comm <- predict(model_new,xx$DATA[comm,])          #无差异的基因
    # summary(pred_comm)
    # comm_error <- 1-sum(pred_comm)/nrow(xx$DATA[comm,])
    # comm_error
    # # 
    #所有有差异的基因
    pred_test <- predict(model_new, xx$DATA[diff,])
    # summary(pred_test)
    # 
    # test_error[i] <- sum(pred_test)/nrow(x_test)
    # print(test_error[i])
    
    
    #pred_all <- predict(model_new, xx$DATA)                  #总的数据
    # summary(pred_all)
    # 
    # fdr_svm_sig[i] <- sum(pred_comm == FALSE)/sum(pred_all == FALSE)
    
    test_error[i] <- (sum(pred_comm == FALSE) + sum(pred_test == TRUE))/nrow(xx$DATA)
    
    
}

index <- which(seq(10000)%%50 == 0)


par(mfrow=c(1,1))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))

plot(train_error[c(1,5,10,20,30,40,index)],ylim = c(0,0.5),col="blue",
     xlab = expression(paste("Index of ",symbol(n))), 
     ylab = "Error Rate",pch = 1)
points(test_error[c(1,5,10,20,30,40,index)],col="red",pch = 1)
abline(h=0.1,col="black",lty=2,lwd=2)

legend("topright",c("Test Errors","Training Errors"),
       col = c("red","blue"),pch = 1)










#################################################################
#----MEB AUC-------#
check <- numeric(nrow(xx$DATA))
for(m in 1:nrow(xx$DATA)){
    check[m] <- decision_function(xx$DATA[m,], model=model,gamma=gamma[gamma_num_new])-model$rho
}
sum(check>0)           #no.TRUE   non-DE genes
summary(predict(model,xx$DATA))


check_ord_MEB <- order(check)
category <- c(rep(1,length(dii)),rep(0,length(comm)),rep(1,nrow(xx$DATA)-length(ci)))   #DE genes equal to 1
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

roc(category, check)

#################################################
#ROC curve
par(mfrow=c(1,1))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))

#MEB
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

#plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE, asp=NA, ann = F)
plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE,ylim = c(0,1))

par(pty="s")