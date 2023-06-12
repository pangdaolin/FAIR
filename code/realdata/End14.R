library(plyr)
library(caret)
library(textir)
library(randomForest)
library(glmnet)
library(MGLM)
library(mgcv)
library(readr)
library(openxlsx)
source("FAIR.R")
lc_study_otu_table <- read_tsv("/data/lc_study_otu_table.tsv")
lc_study_mapping_file <- read_tsv("/data/lc_study_mapping_file.tsv")
data.id <- read.xlsx("/data/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("/data/pbio.2003862.s040.xlsx")
data1.id <- data.id[which(data.id$Compartment=="Endosphere"&data.id$Season=="2014"),1]
imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Endosphere"),1]

OTU.id <- lc_study_otu_table$OTUID
all.data.id <- colnames(lc_study_otu_table)
X.train1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data1.id,all.data.id)[-which(match(data1.id,all.data.id,nomatch = 0)==0)]]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
X.train1 <- t(as(X.train1,"matrix"))
Y.train1 <- lc_study_mapping_file$Age[match(data1.id,all.data.id)[-which(match(data1.id,all.data.id,nomatch = 0)==0)]]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
X.train1x <- X.train1[,c(1:(22-1),(22+1):(ncol(X.train1)),22)]

X.data <- as.matrix(X.train1x)
Y.data <- Y.train1 

rmse11 <- numeric()
rmse12 <- numeric()
rmse21 <- numeric()
rmse22 <- numeric()
rmse3 <- numeric()
rmse4 <- numeric()
rmse5 <- numeric()
rmse61 <- numeric()
rmse62 <- numeric()


set.seed(123456)
for (i in 1:100){
  
  train.id <- sample(datasize,round(datasize/3*2))
  X.train <- X.data[train.id,]  
  X.test <- X.data[-train.id,]
  Y.train <- Y.data[train.id]
  Y.test <- Y.data[-train.id]
  
  X.test=X.test[,which(colSums(X.train)!=0)]
  X.train=X.train[,which(colSums(X.train)!=0)]
  p=ncol(X.train)
  m <- rowSums(X.train)
  m1 <- rowSums(X.test)
  
  ###############################################
  fit0 <- FAIR(X.train, scale(Y.train), penalty=NULL, n.factors=1, maxit = 500)
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
  
  model1 <- FAIR(X.train1, scale(Y.train), penalty=NULL, n.factors=1, maxit = 500)
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
  
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM1 <- gam(Y~s(X1)+s(X2),data=data.frame(res1))
  pre11 <- predict(GAM1,test.object)
  GAM2 <- gam(Y~s(X1,X2),data=data.frame(res1))
  pre12 <- predict(GAM2,test.object)
  
  rmse11 <- c(rmse11,sqrt(sum((pre11-Y.test)^2)/length(Y.test)))
  rmse12 <- c(rmse12,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  
  ##################
  model3 <- mnlm(cl=NULL, Y.train, X.train,free=1)#, cv=T, nfold=10
  B <- coef(model3)#, select="1se"
  phi3 <- as.matrix(B[-1,,drop=FALSE])

  z31_train <- as.vector(X.train%*%t(phi3))/m
  z31_test <- as.vector(X.test%*%t(phi3))/m1#ZZ1
  
  res2 = data.frame(Y=Y.train,X1=(z31_train))
  test.object1 <- data.frame(X1=z31_test)
  
  GAM2 <- gam(Y~s(X1),data=data.frame(res2))#+s(m)
  pre3 <- predict(GAM2,test.object1)
  rmse3 <- c(rmse3,sqrt(sum((pre3-Y.test)^2)/length(Y.test)))
  
  ###################
  rfdata.train <- data.frame(Y=Y.train,X=X.train/m)
  model5 <- randomForest(Y~.,data = rfdata.train)   #建模，ntree=j 指的树数,ntree = 10000
  rfdata.test <- data.frame(X=X.test/m1)
  pre5 <- predict(model5,rfdata.test)  
  
  rmse5 <- c(rmse5,sqrt(sum((pre5-Y.test)^2)/length(Y.test)))
  
  
  ############################
  model6 <- cv.glmnet(x=X.train/m,y = Y.train,family = "gaussian")
  pre61 <- predict(model6,X.test/m1, s = "lambda.min")   #预测
  pre62 <- predict(model6,X.test/m1, s = "lambda.1se")   #预测
  
  rmse61 <- c(rmse61,sqrt(sum((pre61-Y.test)^2)/length(Y.test)))
  rmse62 <- c(rmse62,sqrt(sum((pre62-Y.test)^2)/length(Y.test)))
  
}



boxdata21 <- data.frame(err=c(rmse12,rmse3,rmse5,rmse62),
                        method=c(rep("FAIR-GAM",100),rep("MNIR-GAM",100),rep("randomForest",100),
                                 rep("glmnet",100)),
                        data=rep("Endosphere",400))



