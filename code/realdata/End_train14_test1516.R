library(plyr)
library(caret)
library(textir)
library(randomForest)
library(glmnet)
library(MGLM)
library(mgcv)
library(readr)
library(openxlsx)

source("/code/FAIR.R")

lc_study_otu_table <- read_tsv("/data/lc_study_otu_table.tsv")
lc_study_mapping_file <- read_tsv("/data/lc_study_mapping_file.tsv")
data.id <- read.xlsx("/data/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("/data/pbio.2003862.s040.xlsx")


data.train1.id <- data.id[which(data.id$Compartment=="Endosphere"&data.id$Season=="2014"),1]
data.test1.id <- data.id[which(data.id$Compartment=="Endosphere"&data.id$Season=="2015"),1]
# data.test1.id <- data.id[which(data.id$Compartment=="Endosphere"&data.id$Season=="2016"),1]
imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Endosphere"),1]
OTU.id <- lc_study_otu_table$OTUID
all.data.id <- colnames(lc_study_otu_table)
X.train1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data.train1.id,all.data.id)[-which(match(data.train1.id,all.data.id,nomatch = 0)==0)]]
X.train1 <- t(as(X.train1,"matrix"))
X.test1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data.test1.id,all.data.id)[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]]
X.test1 <- t(as(X.test1,"matrix"))
Y.train1 <- lc_study_mapping_file$Age[match(data.train1.id,all.data.id)[-which(match(data.train1.id,all.data.id,nomatch = 0)==0)]]
Y.test1 <- lc_study_mapping_file$Age[match(data.test1.id,all.data.id)[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]]

X.train1x <- X.train1[,c(1:(22-1),(22+1):(ncol(X.train1)),22)]
X.test1x <- X.test1[,c(1:(22-1),(22+1):(ncol(X.test1)),22)]


X.train <- X.train1x
X.test <- X.test1x
Y.train <- Y.train1
Y.test <- Y.test1

datasize <- nrow(X.train)
p=ncol(X.train)
m <- rowSums(X.train)
m1 <- rowSums(X.test)

k2 <- ncol(as.matrix(Y.train))
nn <- nrow(as.matrix(Y.train))

fit0 <- FAIR(X.train, scale(Y.train), penalty=T, n.factors=1, maxit = 500)
phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
Beta0 <- c(fit0[["model.coefs"]][["Beta"]],0)
Q=diag(p)-1/p
QC1=cbind(Q%*%(phi0),Q%*%(Beta0))
newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))

if (newbase1 != (p)&newbase1!=1){
  X.traint <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
}
if(newbase1==1){
  X.traint <- X.train[,c((newbase1+1):(p),newbase1)]
}

model1 <- FAIR(X.traint, scale(Y.train), penalty=T, n.factors=1,ngridpt=100,keep.path=T, maxit = 500)
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
GAM2 <- gam(Y~s(X1,X2),data=data.frame(res1))
pre12 <- predict(GAM2,test.object)
err11 <- sqrt(sum((pre12-Y.test)^2)/length(Y.test))

#######################3
rfdata.train <- data.frame(Y=Y.train,X=X.train/rowSums(X.train))
model5 <- randomForest(Y~.,data = rfdata.train) 
rfdata.test <- data.frame(X=X.test/m1)
pre5 <- predict(model5,rfdata.test)  
err12 <- sqrt(sum((pre5-Y.test)^2)/length(Y.test))

###################
model3 <- mnlm(cl=NULL, Y.train, X.train, free=1)
B <- coef(model3)
phi3 <- as.matrix(B[-1,,drop=FALSE])
z31_train <- as.vector(X.train%*%t(phi3))/m
z31_test <- as.vector(X.test%*%t(phi3))/m1#ZZ1
res2 = data.frame(Y=Y.train,X1=(z31_train))
test.object1 <- data.frame(X1=z31_test)
GAM2 <- gam(Y~s(X1),data=data.frame(res2))#+s(m)
pre3 <- predict(GAM2,test.object1)
err13 <- sqrt(sum((pre3-Y.test)^2)/length(Y.test))

###############
model6 <- cv.glmnet(x=X.train/m,y = Y.train,family = "gaussian")
pre62 <- predict(model6,X.test/m1, s = "lambda.1se") 
err14 <- sqrt(sum((pre62-Y.test)^2)/length(Y.test))

#############      output      ################
errdata <- data.frame(Error=c(err11,err12,err13,err14),Methods=c("FAIR-GAM","randomForest","MNIR-GAM","glmnet"))
errdata

