#biological replicates
rm(list = ls())
library('compcodeR')
source("functions.R")



#------------------------------------------------------------------------------
#生成第一个数据集
n.vars1 <- 15000
samples.pc.comm <- 2
samples.pc1 <- samples.pc.comm
n.diffexp1 <- n.vars1*0.6
fraction.upregulated1 <- 1
effect.size1 <- 8


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



B1


summary(B1@variable.annotations$truelog2foldchanges)

ndiff1 <- sum(B1@variable.annotations$differential.expression == 1)    #第一个数据集中差异表达基因的个数

identical(which(B1@variable.annotations$differential.expression == 1), 
          1:ndiff1)


ndata1 <- nrow(B1@count.matrix)      #第一个数据集的总个数


commDE1 <- B1@count.matrix[ndiff1:ndata1,]    #第一个数据集中非差异表达基因的数据
diffDE1 <- B1@count.matrix[1:ndiff1,]         #第一个数据集中差异表达基因的数据
dat1 <- B1@count.matrix                        #第一个数据集数据



#------------------------------------------------------------------------------
n.vars2 <- 10000
samples.pc2 <- samples.pc.comm
n.diffexp2 <- n.vars2*0.4
fraction.upregulated2 <- 0.5
effect.size2 <- 8


B2 <- generateSyntheticData(dataset = "B2", n.vars = n.vars2, 
                            samples.per.cond = samples.pc2, 
                            n.diffexp = n.diffexp2, 
                            repl.id = 1, seqdepth = 1e6,
                            minfact = 0.7, maxfact = 1.4,
                            relmeans = "auto",
                            dispersions = "auto",
                            fraction.upregulated = fraction.upregulated2, 
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



B2


summary(B2@variable.annotations$truelog2foldchanges)

ndiff2 <- sum(B2@variable.annotations$differential.expression == 1)      #第二个数据集中差异表达基因的个数

identical(which(B2@variable.annotations$differential.expression == 1), 
          1:ndiff2)


ndata2 <- nrow(B2@count.matrix)                #第二个数据集的总数


commDE2 <- B2@count.matrix[ndiff2:ndata2,]
diffDE2 <- B2@count.matrix[1:ndiff2,]
dat2 <- B2@count.matrix



#------------------------------------------------------------------------------
#将上面两个数据集合并













#-----------------------------NIMEB--------------------------------------------
id1 <- sample((ndiff1+1):ndata1, 500, replace = FALSE)
id2 <- sample((ndiff2+1):ndata2, 500, replace = FALSE)
x_train <- rbind(dat1[id1,],dat2[id2,])


library(e1071)
gamma <- seq(0,1e-2,1e-5)
train_error <- numeric(length(gamma))
for (i in 1:length(gamma)) {
    model <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification", 
                 kernel = "radial", gamma = gamma[i],
                 nu = 0.01, tolerance = 0.001, 
                 shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
                 na.action = na.omit)
    
    train_error[i] <- 1 - model$tot.accuracy/100
}


#######################################################
gamma_num_new <- which.min(abs(train_error - 0.1))
train_error[gamma_num_new]


model <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification", 
             kernel = "radial", gamma = gamma[gamma_num_new],
             nu = 0.01, tolerance = 0.001, 
             shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
             na.action = na.omit)

summary(model)
pred_train <- predict(model, x_train)
summary(pred_train) 
train_error_all <- 1 - sum(pred_train)/nrow(x_train)
train_error_all

###########################################
#评估
commDE <- rbind(commDE1, commDE2)
pred_comm <- predict(model,commDE)          #无差异的基因
summary(pred_comm)
comm_error <- 1 - sum(pred_comm)/nrow(commDE)
comm_error




#x_test <- xx$DATA[diff,]                            #所有有差异的基因
diffDE <- rbind(diffDE1, diffDE2)
pred_test <- predict(model, diffDE)
summary(pred_test)
test_error <- sum(pred_test)/nrow(diffDE)
test_error



dat <- rbind(dat1, dat2)
pred_all <- predict(model, dat)                  #总的数据
summary(pred_all)




fdr_svm <- sum(pred_comm == FALSE)/sum(pred_all == FALSE)
fdr_svm





#----NIMEB AUC-------#
library(pROC)
check <- numeric(nrow(dat))
for(i in 1:nrow(dat)){
    check[i] <- decision_function(dat[i,], model=model,gamma=model$gamma)-(model)$rho
}


sum(check>0)           #no.TRUE   non-DE genes
check_ord_MEB <- order(check)



category <- c(B1@variable.annotations$differential.expression, 
              B2@variable.annotations$differential.expression)
roc_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
auc_MEB <- auc(roc_MEB)
auc_MEB

#roc(category, check)












count.matrix <- dat
class <- B1@sample.annotations$condition
gg <- class





#------------------------------HTN---------------------------------------------
hkeep <- c(id1, id2 + nrow(dat1))
nn <- ncol(count.matrix)
fww <- rep(1,nn)
fww[1] = 1
for (t in 2:nn) {
    fww[t] <- SCBNM(count.matrix[,c(1,t)], hkeep, a = 0.05)
}


lambda2 <- rowSums(count.matrix)/sum(count.matrix)
cs <- colSums(count.matrix)
effM3 <- cs*fww
expMeanAdj3 <- outer(lambda2, effM3)

exactPadj3 <-exactTestPoisson(dataMatrix = count.matrix, 
                              group1Ind = which(gg == 1),
                              group2Ind = which(gg == 2),
                              meanMatrix = expMeanAdj3, verbose = TRUE)



fdr_htn <- sum(exactPadj3[c((ndiff1+1):ndata1,(ndata1+ndiff2+1):(ndata1+ndata2))] < 0.05)/sum(exactPadj3 < 0.05)

fdr_htn


roc_HTN <- roc(category, exactPadj3)
auc_HTN <- auc(roc_HTN)
auc_HTN




#------------------------------edgeR--------------------------------------------
library(edgeR)
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(class))
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")


fdr_edgeR <- sum(edgeR.pvalues[c((ndiff1+1):ndata1,(ndata1+ndiff2+1):(ndata1+ndata2))] < 0.05)/sum(edgeR.pvalues < 0.05)
fdr_edgeR


roc_edgeR <- roc(category, edgeR.pvalues)
auc_edgeR <- auc(roc_edgeR)
auc_edgeR



#---------------------------Library size---------------------------------------
pm <- Poisson.model(countMatrix = count.matrix, group1 = which(gg == 1), 
                        group2 = which(gg == 2), calcFactor = FALSE)
# pmAdj <- Poisson.model.new(countMatrix = count.matrix, group1 = which(gg == 1), 
#                            group2 = which(gg == 2), calcFactor = TRUE)


fdr_libsize1 <- sum(pm$stats$pval[c((ndiff1+1):ndata1,(ndata1+ndiff2+1):(ndata1+ndata2))] < 0.05)/sum(pm$stats$pval < 0.05)
fdr_libsize1




# fdr_libsize2 <- length(intersect(which(pmAdj$stats$pval < 0.05),
#                                  (ndiff:ndata)))/sum(pmAdj$stats$pval < 0.05)
# fdr_libsize2



roc_libsize1 <- roc(response = category, predictor = pm$stats$pval)
auc_libsize1 <- auc(roc_libsize1)
auc_libsize1


# roc_libsize2 <- roc(category, pmAdj$stats$pval)
# auc(roc_libsize2)



rownames(count.matrix) <- NULL 
#------------------------------DESeq------------------------------------------
library(DESeq)
DESeq.cds = newCountDataSet(countData = count.matrix,
                            conditions = factor(class))
DESeq.cds = estimateSizeFactors(DESeq.cds)
DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum",
                                method = "pooled", fitType = "local")
DESeq.test = nbinomTest(DESeq.cds, "1", "2")
DESeq.pvalues = DESeq.test$pval
DESeq.adjpvalues = p.adjust(DESeq.pvalues, method = "BH")


fdr_deseq <- sum(DESeq.pvalues[c((ndiff1+1):ndata1,(ndata1+ndiff2+1):(ndata1+ndata2))] < 0.05)/sum(DESeq.pvalues < 0.05)
fdr_deseq



roc_DESeq <- roc(category, DESeq.pvalues)
auc_DESeq <- auc(roc_DESeq)
auc_DESeq










#------------------------------NOIseq------------------------------------------
library(NOISeq)
myfactors = data.frame(condsAB = class)
myData <- NOISeq::readData(data = count.matrix, factors = myfactors)
mynoiseq <- noiseqbio(input = myData, k = 0.5, norm = "tmm", factor = "condsAB")

noiseq_prob <- mynoiseq@results[[1]][,"prob"]


fdr_noiseq <- sum(noiseq_prob[c((ndiff1+1):ndata1,(ndata1+ndiff2+1):(ndata1+ndata2))] > 0.8)/sum(noiseq_prob > 0.8)
fdr_noiseq

roc_NOISeq <- roc(category, noiseq_prob)
auc_NOISeq <- auc(roc_NOISeq)
auc_NOISeq




#-----------------------------------------------------------------------------
fdr_svm
fdr_htn
fdr_edgeR
#calcFactor = FALSE
fdr_libsize1
#calcFactor = TRUE
#fdr_libsize2
fdr_deseq
fdr_noiseq

#------------------------------------------------------------------------------
auc_MEB
auc_libsize1
auc_HTN
auc_edgeR
auc_DESeq
auc_NOISeq




#--------------------------------ROC curve-------------------------------------
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))
#require("RColorBrewer")
#miscolores = brewer.pal(8,"Set3")[1:8]

#MEB
roc_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = TRUE)
auc_MEB <- auc(roc_MEB)
auc_MEB

plot(roc_MEB,col="red",legacy.axes=T,print.auc=F,grid =T)


#edgeR
roc_edgeR <- roc(category, edgeR.pvalues,smooth = TRUE)
auc(roc_edgeR)

plot(roc_edgeR,col="yellow",add=T,legacy.axes=T,print.auc=F)


#Libsize
roc_libsize1 <- roc(category, pm$stats$pval,smooth = TRUE)
auc_libsize1 <- auc(roc_libsize1)

plot(roc_libsize1,col="green",add=T,legacy.axes=T,print.auc=F)


#HTN
roc_HTN <- roc(category, exactPadj3,smooth = TRUE)
auc(roc_HTN)

plot(roc_HTN,col="blue",add=T,legacy.axes=T,print.auc=F)





#DESeq
roc_DESeq <- roc(category, DESeq.pvalues,smooth = TRUE)
auc(roc_DESeq)

plot(roc_DESeq,col="black",add=T,legacy.axes=T,print.auc=F)



#NOISeq
roc_NOISeq <- roc(category, noiseq_prob, smooth = TRUE)
auc(roc_NOISeq)

plot(roc_NOISeq,col="orange",add=T,legacy.axes=F,print.auc=F)


legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESeq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)


