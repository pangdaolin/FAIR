library(ggplot2)
library(ggthemes)
mymin = function(y) ifelse(y <= breaks[1],0, breaks[2])  

fig1data1 <- read.csv("/home/daolinpang/project1/code/review/github/res/nfRhiz14.csv")
fig1data2 <- read.csv("/home/daolinpang/project1/code/review/github/res/nfEnd14.csv")
fig1data <- rbind(fig1data1,fig1data2)

fig1 <- ggplot(fig1data,aes(data,err,fill=method))+
  geom_boxplot()+
  facet_grid(type ~ ., scales = "free")+
  geom_rect(aes(xmin = data , xmax = data ,ymin = mymin(err), ymax = err)) +
  labs(x=NULL, y = "Prediction error")+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14), 
        legend.title = element_text(size=14),strip.text=element_blank())



fig2data1 <- read.csv("/home/daolinpang/project1/code/review/github/res/errRhiz14.csv")
fig2data2 <- read.csv("/home/daolinpang/project1/code/review/github/res/errEnd14.csv")
fig2data <- rbind(fig2data1,fig2data2)

fig2 <- ggplot(fig2data2,aes(data,err,fill=factor(method,level = c("glmnet","MNIR-GAM","FAIR-GAM","randomForest"))))+
  geom_boxplot()+
  labs(x="", y = "Prediction error")+
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_discrete(limits=c("glmnet","MNIR-GAM","FAIR-GAM","randomForest"))+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14), 
        legend.title = element_text(size=14))


