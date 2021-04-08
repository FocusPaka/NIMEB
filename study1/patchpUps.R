rm(list=ls())


library(ggplot2)
library(tidyr)
#View(data)
test <- data
test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
#View(test)


#colnames(test) <- c("Methods","pUp=0.1","pUp=0.3","pUp=0.5")
colnames(test) <- c("Methods","pUp=0.5","pUp=0.7","pUp=0.9")

#View(test)
#扁变长
test_gather <- gather(data=test,
                      key=pUp,
                      value=AUC,
                      -Methods)
#View(test_gather)
p5 <- ggplot(data=test_gather)+
        geom_boxplot(mapping=aes(x=Methods,y=AUC,color=Methods))+
        facet_wrap(~pUp)+
        coord_cartesian(ylim = c(0, 1))+
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))


p1
p2
p3
p4
p5
p6

library(patchwork)
p1+p2+p3+p4+p5+p6+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = "collect",ncol = 2)












#--------------------------------------------------------------------
#View(data)
library(ggplot2)
library(tidyr)
test <- data[,c(1,4,3,2)]

test[which(test[,1]=="MEB"),1] <- "NIMEB"
test[which(test[,1]=="libsize"),1] <- "LibrarySize"
test[which(test[,1]=="LibSize"),1] <- "LibrarySize"
colnames(test) <- c("Methods","pUp=0.5","pUp=0.7","pUp=0.9")
#View(test)
#扁变长
test_gather <- gather(data=test,
                      key=pUp,
                      value=AUC,
                      -Methods)
#View(test_gather)
p6 <- ggplot(data=test_gather)+
    geom_boxplot(mapping=aes(x=Methods,y=AUC,color=Methods))+
    facet_wrap(~pUp)+
    coord_cartesian(ylim = c(0, 1))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))



