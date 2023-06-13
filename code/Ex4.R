library(textir)
library(MGLM)
library(textir)
library(mgcv)
library(LassoSIR)
source("FAIR.R")
source("eval_space.R")

source("/Lin2014/cdmm.R")
dyn.load("/Lin2014/cdmm.so")
n=60; p=50; k=1;
D0 <- numeric()
D1 <- numeric()
D2 <- numeric()
D3 <- numeric()
D4 <- numeric()
D5 <- numeric()

phi.fit1 <- rep(1, p)
phi.fit2 <- rep(1, p)
phi.fit3 <- rep(1, p)
phi.fit4 <- rep(1, p)
phi.fit5 <- rep(1, p)
phi.fit6 <- rep(1, p)

Beta.fit1 <- rep(1, p)
Beta.fit2 <- rep(1, p)
Beta.fit3 <- rep(1, p)
Beta.fit5 <- rep(1, p)


set.seed(12529)
alpha <- rnorm(p-1, mean = 0, sd = 1)
phi <- sample(c(-3:-1,1:3),p-1,replace = T)
ind.phi <- c(1,2,3,4,5,6,7,8,9)
phi[-ind.phi] <- 0
phi[ind.phi] <- c(1,1,1,1,0,0,0,0,0)
bb.seq <- c(seq(from=-3,to=-1,by=1),seq(from=1,to=3,by=1))
Beta <- matrix(sample(bb.seq,k*(p-1),replace = TRUE),k,p-1)
ind.Beta <- c(1,2,3,4,5,6,7,8,9)
Beta[-ind.Beta] <- 0
Beta[ind.Beta] <- c(-1,1,-1,1,1,1,1,1,1)
phi.true <- c(phi,0)
Beta.true <- c(abs(Beta),0)
initbase=p
for (runs in 1:100) {
  
  f <- rnorm(n*k, mean = 0, sd = 1)
  ff <- matrix(f,n,k)
  Y <- rnorm(n, mean = 0, sd = 1)#rbinom(n,1,0.5)
  V = Y
  stru1 <- matrix(alpha, n, p-1, byrow = TRUE) + ff %*% Beta + V %*% t(phi)
  q <- cbind(exp(stru1),1)/(rowSums(exp(stru1))+1)
  
  X <- matrix(0,n,p)
  for (j in 1:n) {
    count <- sample(1:p, (sample(1000:1e4,1)), replace = T, prob = q[j, ])
    for (i in 1:p) {
      X[j,i] <- sum(count==i)
    }
  }
  
  X.train <- X[1:(n/2),];V.train <- V[1:(n/2)];ff.train <- ff[1:(n/2),]
  X.test <- X[(n/2+1):n,];V.test <- V[(n/2+1):n];ff.test <- ff[(n/2+1):n,]
  
  fit10 <- FAIR(X.train, (V.train),base=initbase, penalty=NULL, n.factors=k, maxit = 500)#scale
  phi10 <- c(fit10[["model.coefs"]][["phi"]],0)
  Beta10 <- c(fit10[["model.coefs"]][["Beta"]],0)
  Q=diag(p)-1/p
  # E=rbind(diag(p-1),-1)
  QC1=cbind(Q%*%(phi10),Q%*%(Beta10))
  newbase1=which.min(sqrt(rowSums(QC1^2)))#+(initbase<=which.min(sqrt(rowSums(QC1^2))))
  
  fit1 <- FAIR(X.train, (V.train),base=newbase1, penalty=TRUE, n.factors=k, maxit = 500)#scale
  
  phi1 <- fit1[["model.coefs4"]][["phi"]]
  Beta1 <- fit1[["model.coefs4"]][["Beta"]]
  phi.fit1 <- cbind(phi.fit1, c(phi1,0))
  Beta.fit1 <- cbind(Beta.fit1, c(Beta1,0))
  
  if (newbase1 != (p)&newbase1!=1){
    X.train1 <- X.train[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c(1:(newbase1-1),(newbase1+1):(p),newbase1)]
  }
  if(newbase1==1){
    X.train1 <- X.train[,c((newbase1+1):(p),newbase1)]
    X.test1 <- X.test[,c((newbase1+1):(p),newbase1)]
  }
  
  if (initbase != (p)&initbase!=1){
    X.train0 <- X.train[,c(1:(initbase-1),(initbase+1):(p),initbase)]
  }
  if(initbase==1){
    X.train0 <- X.train[,c((initbase+1):(p),initbase)]
  }
  if(initbase==p){
    X.train0 <- X.train
  }
  
  fits2 <- mnlm(cl=NULL, V.train, X.train0, nlambda=10000, lambda.start=Inf, lambda.min.ratio=0.01)
  B <- coef(fits2, k=log(ncol(X.train)), corrected=FALSE)
  # fits3 <- mnlm(cl=NULL, V.train, X.train, cv=T, nfold=10)
  # B <- coef(fits3, select="min")
  # B2 <- coef(fits3, select="1se")
  phi2 <- as.vector(B[-1,,drop=FALSE])
  phi.fit2 <- cbind(phi.fit2, phi2)
  
  
  
  ##############
  V.train.new1 <- cbind(ff.train,V.train)
  fits30 <- MGLMreg(formula=X.train0~V.train.new1, dist="MN")
  coef.fits3 <- fits30@coefficients
  coef3 <- cbind(coef.fits3[2:3,],0)
  QC3=Q%*%t(coef3)
  newbase3=which.min(sqrt(rowSums(QC3^2)))#+(initbase<=which.min(sqrt(rowSums(QC3^2))))
  
  if (newbase3 != (p)&newbase3!=1){
    X.train3 <- X.train[,c(1:(newbase3-1),(newbase3+1):(p),newbase3)]
    X.test3 <- X.test[,c(1:(newbase3-1),(newbase3+1):(p),newbase3)]
  }
  if(newbase3==1){
    X.train3 <- X.train[,c((newbase3+1):(p),newbase3)]
    X.test3 <- X.test[,c((newbase3+1):(p),newbase3)]
  }
  
  fit3 <- LFIRtune(formula=X.train3~1+V.train.new1, dist="MN", k=k, penalty="group",va.sigma=0, ngridpt=100, penidx=c(FALSE,rep(TRUE,ncol(V.train.new1))))
  goldfit1 <- fit3@select@coefficients
  phi3 <- goldfit1[3,]
  Beta3 <- goldfit1[2,]
  phi.fit3 <- cbind(phi.fit3, c(phi3,0))
  Beta.fit3 <- cbind(Beta.fit3, c(Beta3,0))
  
  
  X.train.zr <- X.train0
  X.train.zr[X.train0 == 0] <- 0.5
  X.train.zr <- X.train.zr/rowSums(X.train.zr)
  logX.train <- log(X.train.zr)
  fit4 <- gic.cdmm(y=V.train, x=logX.train)
  phi4 <- fit4[["bet"]]
  phi.fit4 <- cbind(phi.fit4, phi4)
  
  fit5 <- LassoSIR(X=logX.train,Y=V.train,H=10,no.dim = 2)#choosing.d="manual"
  phi5 <- fit5[["beta"]]
  phi5 <- rowSums(abs(phi5))
  phi.fit5 <- cbind(phi.fit5, phi5)
  
  X.train.a <- X.train0/rowSums(X.train0)
  fit6 <- LassoSIR(X=X.train.a,Y=V.train,H=10,no.dim = 2)#choosing.d="manual"
  phi6 <- fit6[["beta"]]
  phi6 <- rowSums(abs(phi6))
  phi.fit6 <- cbind(phi.fit6, phi6)
  
}


phi.fit1 <- phi.fit1[,-1]#FAIR
phi.fit2 <- phi.fit2[,-1]#MNIR
phi.fit3 <- phi.fit3[,-1]#Oracle
# phi.fit3 <- phi.fit3-matrix(phi.fit3[p,],nrow = p, ncol = 100, byrow = T)
phi.fit4 <- phi.fit4[,-1]#Lin2014
phi.fit5 <- phi.fit5[,-1]#SIR1
phi.fit6 <- phi.fit6[,-1]#SIR2

Beta.fit1 <- (Beta.fit1[,-1])
Beta.fit3 <- (Beta.fit3[,-1])


#####################
fun1 <- function(x){
  which(abs(x)>1e-5)
}
fun2 <- function(x,S.true){
  sum(x %in% S.true)
}
fun3 <- function(x,S.true){
  sum(!x %in% S.true)
}
fun4 <- function(x,S.true){
  all(S.true %in% x)
}
fun5 <- function(x,S.true){
  identical(S.true, x)
}


PB.true <- c(1,2,3,4,5,6,7,8,9)#
# PB.true <- c(S.true,B.true)

PB1.fit <- abs(Beta.fit1)+abs(phi.fit1)
PB2.fit <- abs(phi.fit2)
PB3.fit <- abs(Beta.fit3)+abs(phi.fit3)
PB4.fit <- abs(phi.fit4)
PB5.fit <- abs(phi.fit5)
PB6.fit <- abs(phi.fit6)


PB1.hat <- apply(PB1.fit,2,fun1)
PB2.hat <- apply(PB2.fit,2,fun1)
PB3.hat <- apply(PB3.fit,2,fun1)
PB4.hat <- apply(PB4.fit,2,fun1)
PB5.hat <- apply(PB5.fit,2,fun1)
PB6.hat <- apply(PB6.fit,2,fun1)

TP12 <- lapply(PB1.hat,fun2,S.true=PB.true)
TP22 <- lapply(PB2.hat,fun2,S.true=PB.true)
TP32 <- lapply(PB3.hat,fun2,S.true=PB.true)
TP42 <- lapply(PB4.hat,fun2,S.true=PB.true)
TP52 <- lapply(PB5.hat,fun2,S.true=PB.true)
TP62 <- lapply(PB6.hat,fun2,S.true=PB.true)

FP12 <- lapply(PB1.hat,fun3,S.true=PB.true)
FP22 <- lapply(PB2.hat,fun3,S.true=PB.true)
FP32 <- lapply(PB3.hat,fun3,S.true=PB.true)
FP42 <- lapply(PB4.hat,fun3,S.true=PB.true)
FP52 <- lapply(PB5.hat,fun3,S.true=PB.true)
FP62 <- lapply(PB6.hat,fun3,S.true=PB.true)

# fun5 <- function(x,S.true){
#   identical((S.true), as.numeric(x))
# }


TP12.me <- mean(unlist(TP12))/9
TP12.sd <- sd(unlist(TP12))/9
TP22.me <- mean(unlist(TP22))/9
TP22.sd <- sd(unlist(TP22))/9
TP32.me <- mean(unlist(TP32))/9
TP32.sd <- sd(unlist(TP32))/9
TP42.me <- mean(unlist(TP42))/9
TP42.sd <- sd(unlist(TP42))/9
TP52.me <- mean(unlist(TP52))/9
TP52.sd <- sd(unlist(TP52))/9
TP62.me <- mean(unlist(TP62))/9
TP62.sd <- sd(unlist(TP62))/9


FP12.me <- mean(unlist(FP12))/(p-9)
FP12.sd <- sd(unlist(FP12))/(p-9)
FP22.me <- mean(unlist(FP22))/(p-9)
FP22.sd <- sd(unlist(FP22))/(p-9)
FP32.me <- mean(unlist(FP32))/(p-9)
FP32.sd <- sd(unlist(FP32))/(p-9)
FP42.me <- mean(unlist(FP42))/(p-9)
FP42.sd <- sd(unlist(FP42))/(p-9)
FP52.me <- mean(unlist(FP52))/(p-9)
FP52.sd <- sd(unlist(FP52))/(p-9)
FP62.me <- mean(unlist(FP62))/(p-9)
FP62.sd <- sd(unlist(FP62))/(p-9)

table3data <- data.frame(data=rep(paste("Ex.4 n =",n,"p =",p),3),
                         method=c("FAIR","MNIR","MNIR-f","CLASSO","LASSO-SIR_1","LASSO-SIR_2"),
                         TPR.mean=c(TP12.me,TP22.me,TP32.me,TP42.me,TP52.me,TP62.me),
                         TPR.sd=c(TP12.sd,TP22.sd,TP32.sd,TP42.sd,TP52.sd,TP62.sd),
                         FPR.mean=c(FP12.me,FP22.me,FP32.me,FP42.me,FP52.me,FP62.me),
                         FPR.sd=c(FP12.sd,FP22.sd,FP32.sd,FP42.sd,FP52.sd,FP62.sd))

table3data
