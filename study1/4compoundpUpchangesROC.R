rm(list = ls())


#ROC curve
par(mfrow=c(1,3))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))


#MEB
plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE,lwd=3)

par(pty="s")


#HTN
plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)

#edgeR
plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)

#Libsize
plot(roc_temp_LibSize,col="green",add=T,legacy.axes=T)

#DESeq
plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T,print.auc=F)

title(main = "pUp=0.9")

legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESeq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)

