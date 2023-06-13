library(plyr)
library(caret)
library(textir)
library(randomForest)
library(glmnet)
library(MGLM)
library(mgcv)
library(readr)
library(openxlsx)
lc_study_otu_table <- read_tsv("data/lc_study_otu_table.tsv")
lc_study_mapping_file <- read_tsv("data/lc_study_mapping_file.tsv")
data.id <- read.xlsx("data/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("data/pbio.2003862.s040.xlsx")

# data1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$Site=="Arbuckle"&data.id$Season=="2014"),1]
# imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Rhizosphere"),1]

data1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$Season=="2014"),1]
# data1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$type=="Test"),1]
imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Rhizosphere"),1]

OTU.id <- lc_study_otu_table$OTUID
all.data.id <- colnames(lc_study_otu_table)
X.train1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data1.id,all.data.id)]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
X.train1 <- t(as(X.train1,"matrix"))
# X.test1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data.test1.id,all.data.id)[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]]
# X.test1 <- t(as(X.test1,"matrix"))
Y.train1 <- lc_study_mapping_file$Age[match(data1.id,all.data.id)]#[-which(match(data1.id,all.data.id,nomatch = 0)==0)]
# Y.test1 <- lc_study_mapping_file$Age[match(data.test1.id,all.data.id)[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]]
which(rowSums(X.train1)==0)
which(colSums(X.train1==0)==0)
X.train1x <- X.train1[,c(1:(22-1),(22+1):(ncol(X.train1)),22)]

X.data <- as.matrix(X.train1x)
Y.data <- Y.train1 
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
ck=10
datasize <- nrow(X.data)
cvlist <- CVgroup(k = ck,datasize = datasize,seed = 9708)

rmse11 <- numeric()
rmse12 <- numeric()
rmse21 <- numeric()
rmse22 <- numeric()
rmse3 <- numeric()
rmse4 <- numeric()
rmse5 <- numeric()
rmse61 <- numeric()
rmse62 <- numeric()

phi.fit1 <- rep(1, p)
phi.fit2 <- rep(1, p)
phi.fit3 <- rep(1, p)
phi.fit4 <- rep(1, p)
phi.fit7 <- rep(1, p)

Beta.fit1 <- rep(1, p)
Beta.fit2 <- rep(1, p)
Beta.fit3 <- rep(1, p)
Beta.fit5 <- rep(1, p)

set.seed(123456)
for (i in 1:100){
  # i=1
  # X.train1x <- X.data[-cvlist[[i]],]  #刚才通过cvgroup生成的函数
  # X.test1x <- X.data[cvlist[[i]],]
  # Y.train1 <- Y.data[-cvlist[[i]]]
  # Y.test1 <- Y.data[cvlist[[i]]]
  
  train.id <- sample(datasize,round(datasize/3*2))
  X.train <- X.data[train.id,]  #刚才通过cvgroup生成的函数
  X.test <- X.data[-train.id,]
  Y.train <- Y.data[train.id]
  Y.test <- Y.data[-train.id]
  
  X.test=X.test[,which(colSums(X.train)!=0)]
  X.train=X.train[,which(colSums(X.train)!=0)]
  p=ncol(X.train)
  # colSums(X.train1x)
  
  m <- rowSums(X.train)
  m1 <- rowSums(X.test)
  
  ###############################################
  fit0 <- LMNIRpen1(X.train, scale(Y.train), penalty=NULL, n.factors=1, maxit = 500)
  phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
  Beta0 <- c(fit0[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  # E=rbind(diag(p-1),-1)
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
  
  model1 <- LMNIRpen1(X.train1, scale(Y.train), penalty=T, n.factors=1, maxit = 500)
  phi1 <- model1[["model.coefs"]][["phi"]]
  Beta1 <- model1[["model.coefs"]][["Beta"]]
  phi2 <- model1[["model.coefs4"]][["phi"]]
  Beta2 <- model1[["model.coefs4"]][["Beta"]]
  if(newbase1==p){
    phi1 <- c(phi1,0)
    phi2 <- c(phi2,0)
    Beta1 <- c(Beta1,0)
    Beta2 <- c(Beta2,0)
  }
  if(newbase1 != (p)){
    phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
    phi1[newbase1] <- 0
    phi2[(newbase1+1):p] <- phi2[(newbase1):(p-1)]
    phi2[newbase1] <- 0
    Beta1[(newbase1+1):p] <- Beta1[(newbase1):(p-1)]
    Beta1[newbase1] <- 0
    Beta2[(newbase1+1):p] <- Beta2[(newbase1):(p-1)]
    Beta2[newbase1] <- 0
  }
  
  phi.fit1 <- cbind(phi.fit1, phi1)
  phi.fit2 <- cbind(phi.fit2, phi2)
  Beta.fit1 <- cbind(Beta.fit1, Beta1)
  Beta.fit2 <- cbind(Beta.fit2, Beta2)
  
  z11_train <- X.train%*%(phi1)/m
  z12_train <- X.train%*%(Beta1)/m
  z21_train <- X.train%*%(phi2)/m
  z22_train <- X.train%*%(Beta2)/m
  # z1 <- log2(healthy_impdata_train[,-ncol(healthy_impdata_train)]+1)%*%t(phi1)
  # z2 <- log2(healthy_impdata_train[,-ncol(healthy_impdata_train)]+1)%*%t(Beta1)
  # 
  # plot(scale(z2),(scale(Y.train1)))
  # plot(scale(z1),scale(Y.train1))
  # plot(scale(z2),Y.train1)
  # plot(exp(z1),(Y.train1))
  # plot(exp(z2),(Y.train1))
  z11_test <- X.test%*%(phi1)/m1
  z12_test <- X.test%*%(Beta1)/m1
  z21_test <- X.test%*%(phi2)/m1
  z22_test <- X.test%*%(Beta2)/m1
  
  res1 = data.frame(Y=Y.train,X1=(z11_train),X2=(z12_train))#,m=scale(m)
  # res21 = data.frame(Y=Y.test,X1=(z11_test),X2=(z12_test))
  train.object <- data.frame(X1=z11_train,X2=z12_train)
  test.object <- data.frame(X1=z11_test,X2=z12_test)
  
  GAM1 <- gam(Y~s(X1)+s(X2),data=data.frame(res1))
  # GAM1 <- gam(Y~s(X1)+s(X2)+s(X3),data=data.frame(res1))
  # pre102_train <- predict(GAM1,train.object)
  # plot(Y.train1,pre12)
  pre11 <- predict(GAM1,test.object)
  
  GAM2 <- gam(Y~s(X1,X2),data=data.frame(res1))
  # GAM1 <- gam(Y~s(X1)+s(X2)+s(X3),data=data.frame(res1))
  # pre103_train <- predict(GAM2,train.object)
  # plot(Y.train1,pre12)
  pre12 <- predict(GAM2,test.object)
  # GAM_pre1 <- lm(pre101_test~Y.test1)
  # GAM_pre2 <- lm(pre102_test~Y.test1)
  # GAM_pre3 <- lm(pre103_test~Y.test1)
  # GAM_pre4 <- lm(pre104_test~Y.test1)
  # GAM_pre5 <- lm(pre105_test~Y.test1)
  
  rmse11 <- c(rmse11,sqrt(sum((pre11-Y.test)^2)/length(Y.test)))
  rmse12 <- c(rmse12,sqrt(sum((pre12-Y.test)^2)/length(Y.test)))
  
  
  res1 = data.frame(Y=Y.train,X1=(z21_train),X2=(z22_train))#,m=scale(m)
  # res21 = data.frame(Y=Y.test,X1=(z14_test),X2=(z24_test))
  train.object <- data.frame(X1=z21_train,X2=z22_train)
  test.object <- data.frame(X1=z21_test,X2=z22_test)
  
  GAM1 <- gam(Y~s(X1)+s(X2),data=data.frame(res1))
  # GAM1 <- gam(Y~s(X1)+s(X2)+s(X3),data=data.frame(res1))
  # pre102_train <- predict(GAM1,train.object)
  # plot(Y.train1,pre12)
  pre21 <- predict(GAM1,test.object)
  
  GAM2 <- gam(Y~s(X1,X2),data=data.frame(res1))
  # GAM1 <- gam(Y~s(X1)+s(X2)+s(X3),data=data.frame(res1))
  # pre21_train <- predict(GAM2,train.object)
  # plot(Y.train1,pre12)
  pre22 <- predict(GAM2,test.object)
  
  # GAM_pre1 <- lm(pre101_test~Y.test1)
  # GAM_pre2 <- lm(pre102_test~Y.test1)
  # GAM_pre3 <- lm(pre103_test~Y.test1)
  # GAM_pre4 <- lm(pre104_test~Y.test1)
  # GAM_pre5 <- lm(pre105_test~Y.test1)
  
  rmse21 <- c(rmse21,sqrt(sum((pre21-Y.test)^2)/length(Y.test)))
  rmse22 <- c(rmse22,sqrt(sum((pre22-Y.test)^2)/length(Y.test)))
  
  ##################
  model3 <- mnlm(cl=NULL, Y.train, X.train, free=1)
  # B <- coef(fits3, select="1se")
  B <- coef(model3)
  phi3 <- as.matrix(B[-1,,drop=FALSE])
  phi.fit3 <- cbind(phi.fit3, phi3)
  
  z31_train <- as.vector(X.train%*%t(phi3))/m
  z31_test <- as.vector(X.test%*%t(phi3))/m1#ZZ1
  
  res2 = data.frame(Y=Y.train,X1=(z31_train))
  test.object1 <- data.frame(X1=z31_test)
  
  GAM2 <- gam(Y~s(X1),data=data.frame(res2))#+s(m)
  # plot(Y.train1,pre22)
  pre3 <- predict(GAM2,test.object1)
  # GAM_pre1 <- lm(pre611_test~Y.test1)
  # GAM_pre2 <- lm(pre612_test~Y.test1)
  # GAM_pre3 <- lm(pre613_test~Y.test1)
  # GAM_pre4 <- lm(pre614_test~Y.test1)
  rmse3 <- c(rmse3,sqrt(sum((pre3-Y.test)^2)/length(Y.test)))
  
  
  
  model4 <- mnlm(cl=NULL, (Y.train), X.train, nlambda=10000, lambda.start=Inf, lambda.min.ratio=0.01)
  # fits31 <- mnlm(cl=NULL, V.train, X.train, nlambda=10000, lambda.start=Inf, lambda.min.ratio=0.01)
  B1 <- coef(model4, k=log(ncol(X.train)), corrected=FALSE)
  # B1 <- coef(model41)
  phi4 <- as.matrix(B1[-1,,drop=FALSE])
  phi.fit4 <- cbind(phi.fit4, phi4)
  
  z41_train <- as.vector(X.train%*%t(phi4))/m
  z41_test <- as.vector(X.test%*%t(phi4))/m1
  
  res12 = data.frame(Y=Y.train,X1=(z41_train))#scale
  test.object1 <- data.frame(X1=z41_test)
  
  GAM2 <- gam(Y~s(X1),data=data.frame(res12))#+s(m)
  # plot(Y.train1,pre22)
  pre4 <- predict(GAM2,test.object1)
  # GAM_pre1 <- lm(pre611_test~Y.test1)
  # GAM_pre2 <- lm(pre612_test~Y.test1)
  # GAM_pre3 <- lm(pre613_test~Y.test1)
  # GAM_pre4 <- lm(pre614_test~Y.test1)
  rmse4 <- c(rmse4,sqrt(sum((pre4-Y.test)^2)/length(Y.test)))
  
  ###################
  rfdata.train <- data.frame(Y=Y.train,X=X.train/m)
  model5 <- randomForest(Y~.,data = rfdata.train)   #建模，ntree=j 指的树数,ntree = 10000
  rfdata.test <- data.frame(X=X.test/m1)
  pre5 <- predict(model5,rfdata.test)  
  
  # GAM_pre3 <- lm(pre7_test~Y.test1)
  # summary(GAM_pre3)
  # sqrt(sum((pre7_train-Y.train1)^2)/length(Y.train1))
  rmse5 <- c(rmse5,sqrt(sum((pre5-Y.test)^2)/length(Y.test)))
  
  
  # rfdata.train1 <- data.frame(Y=Y.train,X=X.train/m)
  # model51 <- randomForest(Y~.,data = rfdata.train1,ntree = 10000)   #建模，ntree=j 指的树数
  # rfdata.test1 <- data.frame(X=X.test/m1)
  # pre51 <- predict(model51,rfdata.test1)  
  # 
  # # GAM_pre3 <- lm(pre7_test~Y.test1)
  # # summary(GAM_pre3)
  # # sqrt(sum((pre7_train-Y.train1)^2)/length(Y.train1))
  # rmse51 <- c(rmse51,sqrt(sum((pre51-Y.test)^2)/length(Y.test)))
  
  ############################
  model6 <- cv.glmnet(x=X.train/m,y = Y.train,family = "gaussian")
  pre61 <- predict(model6,X.test/m1, s = "lambda.min")   #预测
  # GAM_pre3 <- lm(pre81_test~Y.test1)
  # summary(GAM_pre3)
  # sqrt(sum((pre81_train-Y.train1)^2)/length(Y.train1))
  pre62 <- predict(model6,X.test/m1, s = "lambda.1se")   #预测
  
  rmse61 <- c(rmse61,sqrt(sum((pre61-Y.test)^2)/length(Y.test)))
  rmse62 <- c(rmse62,sqrt(sum((pre62-Y.test)^2)/length(Y.test)))
  
}

mean(rmse11)
mean(rmse12)
mean(rmse21)
mean(rmse22)
mean(rmse3)
mean(rmse4)
mean(rmse5)
mean(rmse61)
mean(rmse62)
mean(rmse7)

sd(rmse11)
sd(rmse12)
sd(rmse21)
sd(rmse22)
sd(rmse3)
sd(rmse4)
sd(rmse5)
sd(rmse61)
sd(rmse62)
sd(rmse7)

boxdata11 <- data.frame(err=c(rmse12,rmse22,rmse3,rmse4,rmse5,rmse62),
                        method=c(rep("FAIR-GAM",100),rep("FAIR (group lasso)-GAM",100),rep("MNIR-GAM",100),rep("MNIR (gamma lasso)-GAM",100),rep("randomForest",100),
                                 rep("glmnet",100)),
                        data=rep("Rhizosphere",600))


