#! /usr/bin/env Rscript
library(plyr)
library(caret)
library(textir)
library(randomForest)
library(glmnet)
library(MGLM)
library(mgcv)
library(openxlsx)
library(readr)

source("/code/FAIR.R")

lc_study_otu_table <- read_tsv("/data/lc_study_otu_table.tsv")# Please extract the compressed file with the same name first
lc_study_mapping_file <- read_tsv("/data/lc_study_mapping_file.tsv")
data.id <- read.xlsx("/data/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("/data/pbio.2003862.s040.xlsx")

data1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$Season=="2014"),1]
imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Rhizosphere"),1]

OTU.id <- lc_study_otu_table$OTUID
all.data.id <- colnames(lc_study_otu_table)
X.train1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data1.id,all.data.id)]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
X.train1 <- t(as(X.train1,"matrix"))
Y.train1 <- lc_study_mapping_file$Age[match(data1.id,all.data.id)]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
X.train1x <- X.train1[,c(1:(22-1),(22+1):(ncol(X.train1)),22)]
X.data <- as.matrix(X.train1x)
Y.data <- Y.train1 


datasize <- nrow(X.data)
p=ncol(X.data)

rmse11 <- numeric()
rmse12 <- numeric()
rmse21 <- numeric()
rmse22 <- numeric()
rmse31 <- numeric()
rmse32 <- numeric()
rmse41 <- numeric()
rmse42 <- numeric()
rmse51 <- numeric()
rmse52 <- numeric()

nf11 <- numeric()
nf12 <- numeric()
nf13 <- numeric()
nf21 <- numeric()
nf22 <- numeric()
nf23 <- numeric()


phi.fit1 <- rep(1, p)
phi.fit2 <- rep(1, p)
phi.fit3 <- rep(1, p)
phi.fit4 <- rep(1, p)
phi.fit7 <- rep(1, p)

Beta.fit1 <- rep(1, p)
Beta.fit2 <- rep(1, p)
Beta.fit3 <- rep(1, p)
Beta.fit5 <- rep(1, p)
# boxplot(rmse401)
#cc
set.seed(123456)

for (runs in 1:100){
  # i=1
  train.id <- sample(datasize,round(datasize/3*2))
  X.train <- X.data[train.id,]  
  X.test <- X.data[-train.id,]
  Y.train <- Y.data[train.id]
  Y.test <- Y.data[-train.id]
  p=ncol(X.train)
  m <- rowSums(X.train)
  m1 <- rowSums(X.test)
  
  k2 <- ncol(as.matrix(Y.train))
  nn <- nrow(as.matrix(Y.train))
  hic1 <- numeric();dic11 <- numeric();dic12 <- numeric()
  hic2 <- numeric();dic21 <- numeric();dic22 <- numeric()
  
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), n.factors=1, maxit = 500)
  i=1
  hic1[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
  pd <- 2*fit0[["lc1"]]-2*fit0[["ELBO1"]]
  dic11[i] <- -2*fit0[["lc1"]]+2*(pd+i*(p-1)+(p-1)*k2)
  dic12[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2+pd)
  
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- c(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  QC1=cbind(Q%*%(phi0),Q%*%(Beta0))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  model1 <- FAIR(X.train1, scale(Y.train), n.factors=1, maxit = 500)
  
  phi1 <- model1[["model.coefs"]][["phi"]]
  Beta1 <- model1[["model.coefs"]][["Beta"]]
  
  if(newbase1==p){
    phi1 <- c(phi1,0)
    Beta1 <- c(Beta1,0)
  }
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    Beta1[(newbase1+1):p] <- Beta1[(newbase1):(p-1)]
    Beta1[newbase1] <- 0
  }
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%(Beta1)/m
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%(Beta1)/m1
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  GAM2 <- gam(Y~s(X1,X2),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  rmse11 <- c(rmse11,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), n.factors=2, maxit = 500)
  i=2
  hic1[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
  pd <- 2*fit0[["lc1"]]-2*fit0[["ELBO1"]]
  dic11[i] <- -2*fit0[["lc1"]]+2*(pd+i*(p-1)+(p-1)*k2)
  dic12[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2+pd)
  
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- cbind(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  QC1=cbind(Q%*%(phi0),Q%*%t(Beta0))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  model1 <- FAIR(X.train1, scale(Y.train), n.factors=2, maxit = 500)
  phi1 <- c(model1[["model.coefs"]][["phi"]],0)
  Beta1 <- cbind(model1[["model.coefs"]][["Beta"]],0)
  
  
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    Beta1[,(newbase1+1):p] <- Beta1[,(newbase1):(p-1)]
    Beta1[,newbase1] <- 0
  }
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%t(Beta1)/m
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%t(Beta1)/m1
  
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM2 <- gam(Y~s(X1,X2.1,X2.2,k=30),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  rmse21 <- c(rmse21,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), n.factors=3, maxit = 500)
  i=3
  hic1[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
  pd <- 2*fit0[["lc1"]]-2*fit0[["ELBO1"]]
  dic11[i] <- -2*fit0[["lc1"]]+2*(pd+i*(p-1)+(p-1)*k2)
  dic12[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2+pd)
  
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- cbind(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  QC1=cbind(Q%*%(phi0),Q%*%t(Beta0))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  model1 <- FAIR(X.train1, scale(Y.train), n.factors=3, maxit = 500)
  phi1 <- c(model1[["model.coefs"]][["phi"]],0)
  Beta1 <- cbind(model1[["model.coefs"]][["Beta"]],0)
  
  
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    Beta1[,(newbase1+1):p] <- Beta1[,(newbase1):(p-1)]
    Beta1[,newbase1] <- 0
  }
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%t(Beta1)/m
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%t(Beta1)/m1
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM2 <- gam(Y~s(X1,X2.1,X2.2,X2.3,k=40),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  rmse31 <- c(rmse31,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), n.factors=4, maxit = 500)
  i=4
  hic1[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
  pd <- 2*fit0[["lc1"]]-2*fit0[["ELBO1"]]
  dic11[i] <- -2*fit0[["lc1"]]+2*(pd+i*(p-1)+(p-1)*k2)
  dic12[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2+pd)
  
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- cbind(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  QC1=cbind(Q%*%(phi0),Q%*%t(Beta0))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  model1 <- FAIR(X.train1, scale(Y.train), n.factors=4, maxit = 500)
  phi1 <- c(model1[["model.coefs"]][["phi"]],0)
  Beta1 <- cbind(model1[["model.coefs"]][["Beta"]],0)
  
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    Beta1[,(newbase1+1):p] <- Beta1[,(newbase1):(p-1)]
    Beta1[,newbase1] <- 0
  }
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%t(Beta1)/m
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%t(Beta1)/m1
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM2 <- gam(Y~s(X1)+s(X2.1)+s(X2.2)+s(X2.3)+s(X2.4),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  rmse41 <- c(rmse41,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), n.factors=5, maxit = 500)
  i=5
  hic1[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
  pd <- 2*fit0[["lc1"]]-2*fit0[["ELBO1"]]
  dic11[i] <- -2*fit0[["lc1"]]+2*(pd+i*(p-1)+(p-1)*k2)
  dic12[i] <- -2*fit0[["lc1"]]+log(nn)*(i*(p-1)+(p-1)*k2+pd)
  
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- cbind(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  QC1=cbind(Q%*%(phi0),Q%*%t(Beta0))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  model1 <- FAIR(X.train1, scale(Y.train), n.factors=5, maxit = 500)
  phi1 <- c(model1[["model.coefs"]][["phi"]],0)
  Beta1 <- cbind(model1[["model.coefs"]][["Beta"]],0)
  
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    Beta1[,(newbase1+1):p] <- Beta1[,(newbase1):(p-1)]
    Beta1[,newbase1] <- 0
  }
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%t(Beta1)/m
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%t(Beta1)/m1
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  # res21 = data.frame(Y=Y.test,X1=(z11_test),X2=(z12_test))
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM2 <- gam(Y~s(X1)+s(X2.1)+s(X2.2)+s(X2.3)+s(X2.4)+s(X2.5),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  rmse51 <- c(rmse51,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  #############################
  
  nf11[runs] <- which.min(hic1)+1
  nf12[runs] <- which.min(dic11)+1
  nf13[runs] <- which.min(dic12)+1
  
  print(paste("iterations=",runs,";err1=",mean(rmse11),";err2=",mean(rmse21),";err3=",mean(rmse31),";err4=",mean(rmse41),";err5=",mean(rmse51)))
}



boxdata11 <- data.frame(err=c(rmse11,rmse21,rmse31,rmse41,rmse51),
                        method=c(rep("d=1",100),rep("d=2",100),rep("d=3",100),
                                 rep("d=4",100),rep("d=5",100)),
                        data=rep("Rhizosphere",500))

write.table(boxdata11,"/code/realdata/nfRhiz14.csv",row.names=FALSE,col.names=TRUE,sep=",")

