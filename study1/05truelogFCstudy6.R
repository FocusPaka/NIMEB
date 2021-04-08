rm(list = ls())
library('compcodeR')


n.vars <- 15000
samples.pc <- 2
pD <- 0.2
n.diffexps <- n.vars * pD
fraction.upregulated <- 0.9
effect.size <- 8
dispersions <- rep(0.1,n.vars)
B_625_625 <- generateSyntheticData(dataset = "B_625_625", 
                                   n.vars = n.vars, 
                                   samples.per.cond = samples.pc, 
                                   n.diffexp = n.diffexps, 
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


truelogfoldchangesTotal <- B_625_625@variable.annotations$truelog2foldchanges



x <- (1:length(truelogfoldchangesTotal))[which(truelogfoldchangesTotal ==0)]
y <- truelogfoldchangesTotal[which(truelogfoldchangesTotal==0)]

x1 <- (1:length(truelogfoldchangesTotal))[which(truelogfoldchangesTotal!=0)]
y1 <- truelogfoldchangesTotal[which(truelogfoldchangesTotal!=0)]

plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("true"," ",log[2],"FC")),main = "Study 6")
points(x1,y1,col="red",pch = 16,ylim = c(-4,4),xlim = c(0,15000))

