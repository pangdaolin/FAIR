library(plyr)
library(caret)
library(textir)
library(randomForest)
library(glmnet)
library(MGLM)
library(mgcv)

library(readr)
library(openxlsx)
cal_2016 <- readRDS("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/cal_2016.rds")
# length(unique(cal_2016[,1]))
cal_ark_data <- readRDS("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/cal_ark_data.rds")
gg_otus_tax <- readRDS("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/gg_otus_tax.rds")
organelle <- readRDS("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/organelle.rds")
lc_study_otu_table <- read_tsv("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/lc_study_otu_table.tsv")
lc_study_mapping_file <- read_tsv("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/lc_study_mapping_file.tsv")
data.id <- read.xlsx("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/pbio.2003862.s004.xlsx")	
imp.otu <- read.xlsx("/Users/pangdaolin/Documents/mypro/2022/mypro/2021/data/plant age/pbio.2003862.s040.xlsx")


data.train1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$Season=="2014"),1]#data.id$Site=="Arbuckle"&
data.test1.id <- data.id[which(data.id$Compartment=="Rhizosphere"&data.id$Season=="2016"),1]#data.id$Site=="Arbuckle"&
imp.otu1 <- imp.otu[which(imp.otu$Compartment=="Rhizosphere"),1]
OTU.id <- lc_study_otu_table$OTUID
all.data.id <- colnames(lc_study_otu_table)
X.train1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data.train1.id,all.data.id)]
X.train1 <- t(as(X.train1,"matrix"))
X.test1 <- lc_study_otu_table[match(imp.otu1,OTU.id),match(data.test1.id,all.data.id)]#[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]
X.test1 <- t(as(X.test1,"matrix"))
Y.train1 <- lc_study_mapping_file$Age[match(data.train1.id,all.data.id)]
Y.test1 <- lc_study_mapping_file$Age[match(data.test1.id,all.data.id)]#[-which(match(data.test1.id,all.data.id,nomatch = 0)==0)]

# which(rowSums(X.train1)==0)
# which(colSums(X.train1==0)==0)
X.train1x <- X.train1[,c(1:(22-1),(22+1):(ncol(X.train1)),22)]
X.test1x <- X.test1[,c(1:(22-1),(22+1):(ncol(X.test1)),22)]

# X.data <- rbind(X.train1x,X.test1x)
# # X.data <- X.data/rowSums(X.data)*round(median(rowSums(X.data)))#相对丰度
# 
# Y.data <- c(Y.train1,Y.test1)

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

fit0 <- LMNIRpen1(X.train, scale(Y.train), penalty=T, n.factors=1, maxit = 500)
phi0 <- c(fit0[["model.coefs"]][["phi"]],0)
Beta0 <- c(fit0[["model.coefs"]][["Beta"]],0)
Q=diag(p)-1/p
# E=rbind(diag(p-1),-1)
QC1=cbind(Q%*%(phi0),Q%*%(Beta0))
newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))

if (newbase1 != (p)&newbase1!=1){
  X.traint <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
}
if(newbase1==1){
  X.traint <- X.train[,c((newbase1+1):(p),newbase1)]
}

model1 <- LMNIRpen1(X.traint, scale(Y.train), penalty=T, n.factors=1,ngridpt=100,keep.path=T, maxit = 500)

phi1 <- model1[["model.coefs"]][["phi"]]
Beta1 <- model1[["model.coefs"]][["Beta"]]
phi2 <- model1[["model.coefs4"]][["phi"]]
Beta2 <- model1[["model.coefs4"]][["Beta"]]



phi_top5 <- model1[["groupsolution"]]@select.list[[46]]@coefficients[2,]
Beta_top5 <- model1[["groupsolution"]]@select.list[[46]]@coefficients[3,]

phi_top10 <- model1[["groupsolution"]]@select.list[[54]]@coefficients[2,]
Beta_top10 <- model1[["groupsolution"]]@select.list[[54]]@coefficients[3,]

phi_top20 <- model1[["groupsolution"]]@select.list[[70]]@coefficients[2,]
Beta_top20 <- model1[["groupsolution"]]@select.list[[70]]@coefficients[3,]

phi_top35 <- model1[["groupsolution"]]@select.list[[79]]@coefficients[2,]
Beta_top35 <- model1[["groupsolution"]]@select.list[[79]]@coefficients[3,]

# phi_top40 <- model1[["groupsolution"]]@select.list[[80]]@coefficients[2,]
# Beta_top40 <- model1[["groupsolution"]]@select.list[[80]]@coefficients[3,]


if(newbase1==p){
  phi1 <- c(phi1,0)
  phi2 <- c(phi2,0)
  phi_top5 <- c(phi_top5,0)
  phi_top10 <- c(phi_top10,0)
  phi_top20 <- c(phi_top20,0)
  phi_top35 <- c(phi_top35,0)
  # phi_top40 <- c(phi_top40,0)
  
  Beta1 <- c(Beta1,0)
  Beta2 <- c(Beta2,0)
  Beta_top5 <- c(Beta_top5,0)
  Beta_top10 <- c(Beta_top10,0)
  Beta_top20 <- c(Beta_top20,0)
  Beta_top35 <- c(Beta_top35,0)
  # Beta_top40 <- c(Beta_top40,0)
}
if(newbase1 != (p)){
  phi1[(newbase1+1):p] <- phi1[(newbase1):(p-1)]
  phi1[newbase1] <- 0
  phi2[(newbase1+1):p] <- phi2[(newbase1):(p-1)]
  phi2[newbase1] <- 0
  phi_top5[(newbase1+1):p] <- phi_top5[(newbase1):(p-1)]
  phi_top5[newbase1] <- 0
  phi_top10[(newbase1+1):p] <- phi_top10[(newbase1):(p-1)]
  phi_top10[newbase1] <- 0
  phi_top20[(newbase1+1):p] <- phi_top20[(newbase1):(p-1)]
  phi_top20[newbase1] <- 0
  phi_top35[(newbase1+1):p] <- phi_top35[(newbase1):(p-1)]
  phi_top35[newbase1] <- 0
  # phi_top40[(newbase1+1):p] <- phi_top40[(newbase1):(p-1)]
  # phi_top40[newbase1] <- 0
  
  
  Beta1[(newbase1+1):p] <- Beta1[(newbase1):(p-1)]
  Beta1[newbase1] <- 0
  Beta2[(newbase1+1):p] <- Beta2[(newbase1):(p-1)]
  Beta2[newbase1] <- 0
  Beta_top5[(newbase1+1):p] <- Beta_top5[(newbase1):(p-1)]
  Beta_top5[newbase1] <- 0
  Beta_top10[(newbase1+1):p] <- Beta_top10[(newbase1):(p-1)]
  Beta_top10[newbase1] <- 0
  Beta_top20[(newbase1+1):p] <- Beta_top20[(newbase1):(p-1)]
  Beta_top20[newbase1] <- 0
  Beta_top35[(newbase1+1):p] <- Beta_top35[(newbase1):(p-1)]
  Beta_top35[newbase1] <- 0
  # Beta_top40[(newbase1+1):p] <- Beta_top40[(newbase1):(p-1)]
  # Beta_top40[newbase1] <- 0
}

z11_train <- X.train%*%(phi1)/m
z12_train <- X.train%*%(Beta1)/m
z21_train <- X.train%*%(phi2)/m
z22_train <- X.train%*%(Beta2)/m

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

sqrt(sum((pre11-Y.test)^2)/length(Y.test))
sqrt(sum((pre12-Y.test)^2)/length(Y.test))


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

sqrt(sum((pre21-Y.test)^2)/length(Y.test))
sqrt(sum((pre22-Y.test)^2)/length(Y.test))




imp0=1:85
PB.fit=abs(phi2)+abs(Beta2)
imp1=which((PB.fit)!=0)
imp1=intersect(imp0,imp1)

PB.fit_top5=abs(phi_top5)+abs(Beta_top5)
imp_top5=which((PB.fit_top5)!=0)
imp_top5=intersect(imp0,imp_top5)

PB.fit_top10=abs(phi_top10)+abs(Beta_top10)
imp_top10=which((PB.fit_top10)!=0)
imp_top10=intersect(imp0,imp_top10)

PB.fit_top20=abs(phi_top20)+abs(Beta_top20)
imp_top20=which((PB.fit_top20)!=0)
imp_top20=intersect(imp0,imp_top20)

PB.fit_top35=abs(phi_top35)+abs(Beta_top35)
imp_top35=which((PB.fit_top35)!=0)
imp_top35=intersect(imp0,imp_top35)

# PB.fit_top40=abs(phi_top35)+abs(Beta_top35)
# imp_top40=which((PB.fit_top35)!=0)
# imp_top40=intersect(imp0,imp_top35)


length(imp1)
length(imp_top5)
length(imp_top10)
length(imp_top20)
length(imp_top35)

# which((PB.fit)!=0)

# rfdata.train <- data.frame(Y=Y.data,X=X.data)
# model5 <- randomForest(Y~.,data = rfdata.train,ntree = 10000)   #建模，ntree=j 指的树数

rfdata.train <- data.frame(Y=Y.train,X=X.train/rowSums(X.train))
model5 <- randomForest(Y~.,data = rfdata.train)   #建模，ntree=j 指的树数,ntree = 10000
rfdata.test <- data.frame(X=X.test/m1)
pre5 <- predict(model5,rfdata.test)  

# GAM_pre3 <- lm(pre7_test~Y.test1)
# summary(GAM_pre3)
# sqrt(sum((pre7_train-Y.train1)^2)/length(Y.train1))
sqrt(sum((pre5-Y.test)^2)/length(Y.test))

inter_top5 <- intersect(imp_top5,order(model5[["importance"]],decreasing=T)[1:5])
inter_top10 <- intersect(imp_top10,order(model5[["importance"]],decreasing=T)[1:10])
inter_top20 <- intersect(imp_top20,order(model5[["importance"]],decreasing=T)[1:20])
inter_top35 <- intersect(imp_top35,order(model5[["importance"]],decreasing=T)[1:35])
inter_top55 <- intersect(imp1,order(model5[["importance"]],decreasing=T)[1:55])

length(inter_top5)
length(inter_top10)
length(inter_top20)
length(inter_top35)
length(inter_top55)


##################
# model3 <- mnlm(cl=NULL, Y.train, X.train, cv=T, nfold=10)
# # B <- coef(fits3, select="min")
# B <- coef(model3, select="1se")
model3 <- mnlm(cl=NULL, Y.train, X.train, free=1)
# B <- coef(fits3, select="min")
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
sqrt(sum((pre3-Y.test)^2)/length(Y.test))



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
sqrt(sum((pre4-Y.test)^2)/length(Y.test))

model6 <- cv.glmnet(x=X.train/m,y = Y.train,family = "gaussian")
pre61 <- predict(model6,X.test/m1, s = "lambda.min")   #预测
# GAM_pre3 <- lm(pre81_test~Y.test1)
# summary(GAM_pre3)
# sqrt(sum((pre81_train-Y.train1)^2)/length(Y.train1))
pre62 <- predict(model6,X.test/m1, s = "lambda.1se")   #预测

sqrt(sum((pre61-Y.test)^2)/length(Y.test))
sqrt(sum((pre62-Y.test)^2)/length(Y.test))

