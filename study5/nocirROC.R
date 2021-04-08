#biological replicates

rm(list = ls())
library('compcodeR')
source("functions.R")


n.vars <- 15000
samples.pc <- 5
n.diffexp <- n.vars*0.6
fraction.upregulated <- 1
effect.size <- 8
dispersions <- rep(0.1,n.vars)


B_625_625 <- generateSyntheticData(dataset = "B_625_625", n.vars = n.vars, 
                                   samples.per.cond = samples.pc, 
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



B_625_625


summary(B_625_625@variable.annotations$truelog2foldchanges)

ndiff <- sum(B_625_625@variable.annotations$differential.expression == 1)

identical(which(B_625_625@variable.annotations$differential.expression == 1), 
          1:ndiff)


ndata <- nrow(B_625_625@count.matrix)


commDE <- B_625_625@count.matrix[ndiff:ndata,]
diffDE <- B_625_625@count.matrix[1:ndiff,]
dat <- B_625_625@count.matrix


#-----------------------------NIMEB--------------------------------------------
id <- sample(ndiff:ndata, 1000, replace = FALSE)
x_train <- B_625_625@count.matrix[id,]


library(e1071)
gamma <- seq(0,1e-3,1e-6)
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
pred_comm <- predict(model,commDE)          #无差异的基因
summary(pred_comm)
comm_error <- 1 - sum(pred_comm)/nrow(commDE)
comm_error




#x_test <- xx$DATA[diff,]                            #所有有差异的基因
pred_test <- predict(model, diffDE)
summary(pred_test)
test_error <- sum(pred_test)/nrow(diffDE)
test_error




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



category <- B_625_625@variable.annotations$differential.expression
roc_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB],smooth = T)
auc_MEB <- auc(roc_MEB)
auc_MEB

#roc(category, check)












count.matrix <- B_625_625@count.matrix
class <- B_625_625@sample.annotations$condition
gg <- class
#------------------------------HTN---------------------------------------------

hkeep <- id
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

exactPadj3 <-exactTestPoisson(dataMatrix = count.matrix, group1Ind = which(gg == 1),
                              group2Ind = which(gg == 2),
                              meanMatrix = expMeanAdj3, verbose = TRUE)



fdr_htn <- length(intersect(which(exactPadj3 < 0.05),
                            (ndiff:ndata)))/sum(exactPadj3 < 0.05)
fdr_htn


roc_HTN <- roc(category, exactPadj3,smooth = T)
auc_HTN <- auc(roc_HTN)





#------------------------------edgeR--------------------------------------------
library(edgeR)
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(class))
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")


fdr_edgeR <- length(intersect(which(edgeR.pvalues < 0.05),
                              (ndiff:ndata)))/sum(edgeR.pvalues < 0.05)
fdr_edgeR


roc_edgeR <- roc(category, edgeR.pvalues,smooth = T)
auc_edgeR <- auc(roc_edgeR)




#---------------------------Library size---------------------------------------
pm <- Poisson.model(countMatrix = count.matrix, group1 = which(gg == 1), 
                        group2 = which(gg == 2), calcFactor = FALSE)
# pmAdj <- Poisson.model.new(countMatrix = count.matrix, group1 = which(gg == 1), 
#                            group2 = which(gg == 2), calcFactor = TRUE)


fdr_libsize1 <- length(intersect(which(pm$stats$pval < 0.05),
                                 (ndiff:ndata)))/sum(pm$stats$pval < 0.05)
fdr_libsize1




# fdr_libsize2 <- length(intersect(which(pmAdj$stats$pval < 0.05),
#                                  (ndiff:ndata)))/sum(pmAdj$stats$pval < 0.05)
# fdr_libsize2



roc_libsize1 <- roc(category, pm$stats$pval,smooth = T)
auc_libsize1 <- auc(roc_libsize1)



# roc_libsize2 <- roc(category, pmAdj$stats$pval)
# auc(roc_libsize2)




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


fdr_deseq <- length(intersect(which(DESeq.pvalues < 0.05),
                              (ndiff:ndata)))/sum(DESeq.pvalues < 0.05)
fdr_deseq



roc_DESeq <- roc(category, DESeq.pvalues,smooth = T)
auc_DESeq <- auc(roc_DESeq)











#------------------------------NOIseq------------------------------------------
library(NOISeq)
myfactors = data.frame(condsAB = class)
myData <- NOISeq::readData(data = count.matrix, factors = myfactors)
mynoiseq <- noiseqbio(input = myData, k = 0.5, norm = "tmm", factor = "condsAB")


noiseq_prob <- mynoiseq@results[[1]][,"prob"]


fdr_noiseq <- length(intersect(which(noiseq_prob > 0.8),
                               (ndiff:ndata)))/sum(noiseq_prob > 0.8)
fdr_noiseq

roc_NOISeq <- roc(category, noiseq_prob,smooth = T)
auc_NOISeq <- auc(roc_NOISeq)
auc_NOISeq




#-----------------------------------------------------------------------------
fdr_svm
fdr_htn
fdr_edgeR
#calcFactor = FALSE
fdr_libsize1
#calcFactor = TRUE
fdr_libsize2
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
roc_libsize <- roc(category, pm$stats$pval,smooth = TRUE)
auc_libsize <- auc(roc_libsize)

plot(roc_libsize,col="green",add=T,legacy.axes=T,print.auc=F)


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




