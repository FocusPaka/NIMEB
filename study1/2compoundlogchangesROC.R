# Fig 1 compound log2foldchanges=1,2,3

rm(list = ls())
library(pROC)


par(mfrow=c(1,3))
oldpar=par(mar=c(3,2.6,1,0.2),mgp=c(1.7,0.5,0))


plot(roc_temp_MEB,col="red",legacy.axes=T,grid=TRUE, ylim = c(0,1))

par(pty="s")

plot(roc_temp_HTN,col="blue",add=T,legacy.axes=T)


plot(roc_temp_edgeR,col="yellow",add=T,legacy.axes=T)


plot(roc_temp_LibSize,col="green",add=T,legacy.axes=T)


plot(roc_temp_DESeq,col="black",add=T,legacy.axes=T)

plot(roc_temp_NOISeq,col="orange",add=T,legacy.axes=T,print.auc=F)

title(main = expression(paste(log[2],"FC=1")))

legend("bottomright",c("NIMEB","HTN","edgeR","Library Size","DESeq","NOISeq"),
       col = c("red","blue","yellow","green","black","orange"),
       lwd=1, cex=0.8)
