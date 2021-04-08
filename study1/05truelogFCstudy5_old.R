rm(list = ls())
library('compcodeR')



B_625_625 <- generateSyntheticData(dataset = "B_625_625", 
                                   n.vars = 15000, 
                                   samples.per.cond = 2, 
                                   n.diffexp = 15000*0.3, 
                                   repl.id = 1, seqdepth = 1e7,
                                   minfact = 0.7, maxfact = 1.4,
                                   relmeans = "auto",
                                   dispersions = "auto",
                                   fraction.upregulated = 0.9, 
                                   between.group.diffdisp = FALSE, 
                                   filter.threshold.total = 1, 
                                   filter.threshold.mediancpm = 0, 
                                   fraction.non.overdispersed = 1, 
                                   random.outlier.high.prob = 0,
                                   random.outlier.low.prob = 0,
                                   single.outlier.high.prob = 0,
                                   single.outlier.low.prob = 0,
                                   effect.size = 4,
                                   output.file = NULL)

B_625_625


summary(B_625_625@variable.annotations$truelog2foldchanges)


truelogfoldchangesTotal <- B_625_625@variable.annotations$truelog2foldchanges



x <- (1:length(truelogfoldchangesTotal))[which(truelogfoldchangesTotal ==0)]
y <- truelogfoldchangesTotal[which(truelogfoldchangesTotal==0)]

x1 <- (1:length(truelogfoldchangesTotal))[which(truelogfoldchangesTotal!=0)]
y1 <- truelogfoldchangesTotal[which(truelogfoldchangesTotal!=0)]

plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("true"," ",log[2],"FC")),main = "Study 5")
points(x1,y1,col="red",pch = 16,ylim = c(-3,3),xlim = c(0,15000))


#main = expression(B[625]^{625})
