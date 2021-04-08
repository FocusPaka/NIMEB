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
pUp <- 1
pDifferential <- 0.6

xx <- generateDataset(commonTags=15000, uniqueTags=c(1000,500), 
                      foldDifference=foldDiff, pUp=pUp, 
                      pDifferential=pDifferential, 
                      empiricalDist=D[,1], 
                      libLimits=c(.9,1.2)*1e6)
xx$trueFactors


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

gamma <- seq(0.0001,0.005,0.0001)
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
gamma_num_new <- which.min(abs(train_error-0.05))
train_error[gamma_num_new]


model <- svm(x_train, y = NULL, scale = FALSE, type = "one-classification", 
             kernel = "radial", gamma = gamma[gamma_num_new],
             nu = 0.01, tolerance = 0.001, 
             shrinking = TRUE, cross = 5, probability = FALSE, fitted = TRUE,
             na.action = na.omit)

summary(model)
pred_train <- predict(model, x_train)
summary(pred_train) 
train_error_all <- 1-sum(pred_train)/nrow(x_train)
train_error_all

###########################################
#评估
pred_comm <- predict(model,xx$DATA[comm,])          #无差异的基因
summary(pred_comm)
comm_error <- 1-sum(pred_comm)/nrow(xx$DATA[comm,])
comm_error

x_test <- xx$DATA[diff,]                            #所有有差异的基因
pred_test <- predict(model, x_test)
summary(pred_test)
test_error <- sum(pred_test)/nrow(x_test)
test_error

pred_new_data <- predict(model, xx$DATA[dii,])      #成倍数的差异基因
summary(pred_new_data)
error <- sum(pred_new_data)/nrow(xx$DATA[dii,])
error


pred_all <- predict(model, xx$DATA)                  #总的数据
summary(pred_all)

pred_all_num <- numeric(length(pred_all))           #non-DE genes
pred_all_num[which(pred_all == FALSE)] <- 1          #DE genes
new_data_c <- cbind(xx$DATA, pred_all_num)

fdr_svm_sig <- sum(pred_comm == FALSE)/sum(pred_all == FALSE)
fdr_svm_sig

#################################################################
#----MEB AUC-------#
check <- numeric(nrow(xx$DATA))
for(m in 1:nrow(xx$DATA)){
    check[m] <- decision_function(xx$DATA[m,], model=model,gamma=gamma[gamma_num_new])-model$rho
}
sum(check>0)           #no.TRUE   non-DE genes
summary(predict(model,xx$DATA))


#oM <- order(sM)
#wM <- oM %in% diff


check_ord_MEB <- order(check)
NIMEBfp <- check_ord_MEB %in% diff

category <- c(rep(1,length(dii)),rep(0,length(comm)),rep(1,nrow(xx$DATA)-length(ci)))   #DE genes equal to 1


#################################################
#ROC curve
par(mfrow=c(1,1))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))

#MEB
roc_temp_MEB <- roc(category[check_ord_MEB], (-check)[check_ord_MEB])
roc_obj_MEB <- auc(roc_temp_MEB)
roc_obj_MEB

#plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE, asp=NA, ann = F)
plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE,ylim = c(0,1))

par(pty="s")
















##################################################
#HTN
hkeep<-id

# ------------------------------------
#(4) calculate Fisher test after SCBN adjusting
# ------------------------------------
fwwS<-SCBNM(xx$DATA,hkeep,a=0.05)
sMm5S <- sage.test(xx$DATA[,1], xx$DATA[,2], n1=sum(xx$DATA[,1])/sqrt(fwwS), 
                   n2=sum(xx$DATA[,2])*sqrt(fwwS))
oMm4S <- order(sMm5S)
wMm4S <- oMm4S %in% diff

roc_temp_HTN <- roc(category, sMm5S)
auc(roc_temp_HTN)

plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)

















#####################################
#edgeR
fW <- calcNormFactors_new(xx$DATA,logratioTrim=0.3, sumTrim=0.05)[2]

#(2) calculate Fisher test after adjusting
# ------------------------------------
sM <- sage.test(xx$DATA[,1], xx$DATA[,2], n1=sum(xx$DATA[,1])/sqrt(fW), 
                n2=sum(xx$DATA[,2])*sqrt(fW))
oM <- order(sM)
wM <- oM %in% diff


#edgeR
roc_temp_edgeR <- roc(category, sM)
auc(roc_temp_edgeR)

plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)




















#####################################
#LibSize
s <- sage.test(xx$DATA[,1], xx$DATA[,2], n1=sum(xx$DATA[,1]), n2=sum(xx$DATA[,2]))
o <- order(s)
w <- o %in% diff

roc_temp_LibSize <- roc(category, s)
auc(roc_temp_LibSize)

plot(roc_temp_LibSize,col="green",add=T,legacy.axes=T)





















#####################################
#没有生物学重复的时候有问题，其pvalue基本上全部为1
#DESeq
colnames(xx$DATA) <- c("sample1","sample2")
type <- factor(c("sample1","sample2"))
database <- xx$DATA
cds <- newCountDataSet(database,type)

#对于没有生物学重复
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method="blind", fitType="local", sharingMode="fit-only"  )
res <- nbinomTest(cds,"sample1","sample2")

ordpvalue <- order(res$pval)
ws <- ordpvalue %in% diff

#DESeq
roc_temp_DESeq <- roc(category, res$pval)
auc(roc_temp_DESeq)

plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)



















#####################################
#NOISeq
myfactors = data.frame(condsAB = c("sample1","sample2"))
myData=NOISeq::readData(data = database, factors = myfactors)
mynoiseq=noiseq(myData,k=0.5,norm="tmm",factor="condsAB",pnr=0.2,nss=5,v=0.02,
                lc=1,replicates="no")

noiseq_prob <- mynoiseq@results[[1]][,"prob"]

ordnoiseq_prob <- order(-noiseq_prob)

nois <- ordnoiseq_prob %in% diff


#NOISeq
roc_temp_NOISeq <- roc(category, noiseq_prob)
auc(roc_temp_NOISeq)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T,print.auc=F)

#title(main = "d=4")

legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESeq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)












###############################################################################
par(mfrow=c(1,1))
#num <- (sum(pred_all == FALSE) + 200)
num <- ndiff
plot(1:num, cumsum(!NIMEBfp[1:num]), ylim = c(0,300),
     xlab="Number of Genes Selected", ylab="Number of False Discoveries", 
     lwd=2,type="l", col = "red")
lines(1:num, cumsum(!w[1:num]),lwd=2,type="l", col= "black" )
lines(1:num, cumsum(!wM[1:num]),lwd=2,type="l", col= "green" )
lines(1:num, cumsum(!wMm4S[1:num]),lwd=2,type="l", col= "yellow")
lines(1:num, cumsum(!ws[1:num]),lwd=2,type="l", col= "orange")
lines(1:num, cumsum(!nois[1:num]),lwd=2,type="l", col= "blue")
legend("topleft", c("NIMEB","Lib size","edgeR","HTN","DESeq","NOISeq"),
       col=c("red","black","green","yellow","orange","blue"),
       lwd=c(2,2,2,2,2,2), cex=.7)
abline(v = sum(pred_all == FALSE), col = 'black', lty = 2)
arrows(x0 = num - 900, y0 = sum(pred_comm == FALSE) + 700, x1 = num-300, 
       y1 = sum(pred_comm == FALSE)+100, length = 0.25, angle = 30,
       code = 2, col = 'blue', lty = par("lty"),
       lwd = 2)
 
points(sum(pred_all == FALSE), sum(pred_comm == FALSE), pch = 3, lwd = 3,col = 'red')












#该方法需要验证一下代码的正确性

num <- ndiff
plot(1:num, cumsum(!w[1:num]), 
     xlab="Number of Genes Selected", ylab="Number of False Discoveries", 
     lwd=2,type="l", col = "red")
#lines(1:num, cumsum(!w[1:num]),lwd=2,type="l", col= "black" )
lines(1:num, cumsum(!wM[1:num]),lwd=2,type="l", col= "green" )
lines(1:num, cumsum(!wMm4S[1:num]),lwd=2,type="l", col= "yellow")
lines(1:num, cumsum(!ws[1:num]),lwd=2,type="l", col= "orange")
lines(1:num, cumsum(!nois[1:num]),lwd=2,type="l", col= "blue")
legend("topleft", c("NIMEB","Lib size","edgeR","HTN","DESeq","NOISeq"),
       col=c("red","black","green","yellow","orange","blue"),
       lwd=c(2,2,2,2,2,2), cex=.7)
abline(v = sum(pred_all == FALSE), col = 'black', lty = 2)
arrows(x0 = num - 900, y0 = sum(pred_comm == FALSE) + 700, x1 = num-300, 
       y1 = sum(pred_comm == FALSE)+100, length = 0.25, angle = 30,
       code = 2, col = 'blue', lty = par("lty"),
       lwd = 2)

points(sum(pred_all == FALSE), sum(pred_comm == FALSE), pch = 3, lwd = 3,col = 'red')


