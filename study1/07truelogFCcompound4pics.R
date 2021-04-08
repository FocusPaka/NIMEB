rm(list = ls())
library(compcodeR)


par(mfrow=c(2,2))
oldpar=par(mar=c(3,2.8,2,0.5),mgp=c(1.7,0.5,0))



#study1/Study 3
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,15500), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 1 & Study 3")
points(x1,y1,col="red",pch = 1,ylim = c(-4,4),xlim = c(0,15500))


legend("bottomright",c("DE genes","non-DE genes"),
       col = c("red","black"),pch = 1)




#study2/Study 4
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,25000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 2 & Study 4")
points(x1,y1,col="red",pch = 1,ylim = c(-4,4),xlim = c(0,25000))





#study5
plot(x,y,pch=1, ylim = c(-5,5),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 5")
points(x1,y1,col="red",pch = 1,ylim = c(-5,5),xlim = c(0,15000))





#study6
plot(x,y,pch=1, ylim = c(-5,5),xlim = c(0,25000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 6")
points(x1,y1,col="red",pch = 1,ylim = c(-5,5),xlim = c(0,25000))



