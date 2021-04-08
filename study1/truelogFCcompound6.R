rm(list = ls())
library(compcodeR)


par(mfrow=c(3,2))
oldpar=par(mar=c(3,2.8,2,0.4),mgp=c(1.7,0.5,0))



#study1
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,15500), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 1")
points(x1,y1,col="red",pch = 1,ylim = c(-4,4),xlim = c(0,15500))







#study2
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,25000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 2")
points(x1,y1,col="red",pch = 1,ylim = c(-4,4),xlim = c(0,25000))


#study3
plot(x,y,pch=16, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 3")
points(x1,y1,col="red",pch = 16,ylim = c(-3,3),xlim = c(0,15000))


#study4
plot(x,y,pch=1, xlab = "Index of genes", ylab = expression(paste("True"," ",log[2],"FC")),
     xlim = c(0,25000), ylim = c(-4,4),main = "Study 4")
points(x1,y1,col="red",pch = 1,ylim = c(-3,3),xlim = c(0,25000))




#study5
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 5")
points(x1,y1,col="red",pch = 1,ylim = c(-3,3),xlim = c(0,15000))


#study6
plot(x,y,pch=1, ylim = c(-4,4),xlim = c(0,15000), xlab = "Index of genes", 
     ylab = expression(paste("True"," ",log[2],"FC")),main = "Study 6")
points(x1,y1,col="red",pch = 1,ylim = c(-4,4),xlim = c(0,15000))

legend("bottomright",c("DE genes","non-DE genes"),
       col = c("red","black"),pch = 1)

