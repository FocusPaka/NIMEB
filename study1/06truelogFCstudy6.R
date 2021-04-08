#biological replicates
rm(list = ls())
library('compcodeR')
source("functions.R")



#------------------------------------------------------------------------------
#生成第一个数据集
n.vars1 <- 15000
samples.pc.comm <- 2
samples.pc1 <- samples.pc.comm
n.diffexp1 <- n.vars1*0.3
fraction.upregulated1 <- 0.9
effect.size1 <- 4
dispersions1 <- rep(0.01,n.vars1)

B1 <- generateSyntheticData(dataset = "B1", n.vars = n.vars1, 
                            samples.per.cond = samples.pc1, 
                            n.diffexp = n.diffexp1, 
                            repl.id = 1, seqdepth = 1e6,
                            minfact = 0.7, maxfact = 1.4,
                            relmeans = "auto",
                            dispersions = dispersions1,
                            fraction.upregulated = fraction.upregulated1, 
                            between.group.diffdisp = FALSE, 
                            filter.threshold.total = 1, 
                            filter.threshold.mediancpm = 0, 
                            fraction.non.overdispersed = 0, 
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
dispersions2 <- rep(0.01,n.vars2)

B2 <- generateSyntheticData(dataset = "B2", n.vars = n.vars2, 
                            samples.per.cond = samples.pc2, 
                            n.diffexp = n.diffexp2, 
                            repl.id = 1, seqdepth = 1e6,
                            minfact = 0.7, maxfact = 1.4,
                            relmeans = "auto",
                            dispersions = dispersions2,
                            fraction.upregulated = fraction.upregulated2, 
                            between.group.diffdisp = FALSE, 
                            filter.threshold.total = 1, 
                            filter.threshold.mediancpm = 0, 
                            fraction.non.overdispersed = 0, 
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



truelogfoldchangesTotal <- c(B1@variable.annotations$truelog2foldchanges,
                             B2@variable.annotations$truelog2foldchanges)
summary(truelogfoldchangesTotal)
hist(truelogfoldchangesTotal)
plot(truelogfoldchangesTotal)



truelogfoldchangesTotalM <- matrix(0, nrow = length(truelogfoldchangesTotal),ncol = 2)

truelogfoldchangesTotalM[,1] <- truelogfoldchangesTotal
truelogfoldchangesTotalM[truelogfoldchangesTotal !=0,2] <- 1

x <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==0)]
y <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==0),1]

x1 <- (1:nrow(truelogfoldchangesTotalM))[which(truelogfoldchangesTotalM[,2]==1)]
y1 <- truelogfoldchangesTotalM[which(truelogfoldchangesTotalM[,2]==1),1]
plot(x,y,pch=1, xlab = "Index of genes", ylab = expression(paste("True"," ",log[2],"FC")),
     xlim = c(0,25000), ylim = c(-4,4),main = "Study 6")
points(x1,y1,col="red",pch = 1,ylim = c(-3,3),xlim = c(0,25000))








