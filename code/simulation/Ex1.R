library(textir)
library(mgcv)
library(MGLM)
library(MASS)
library(CondIndTests)
library(caret)
source("FAIR.R")
source("eval_space.R")

n=60; p=10; k=1;
D0 <- numeric()
D1 <- numeric()
D2 <- numeric()
D3 <- numeric()

phi.fit1 <- rep(1, p)
phi.fit2 <- rep(1, p)
phi.fit3 <- rep(1, p)

Beta.fit1 <- rep(1, p)
Beta.fit3 <- rep(1, p)

eval1 <- c(1,1,1)
eval2 <- c(1,1,1)
eval3 <- c(1,1,1)



set.seed(927017)
alpha <- rep(0, p-1)
phi <- c(1,1,1,1,0,0,0,0,0)
Beta <- c(-1,1,-1,1,1,1,1,1,1)
system.time(for (runs in 1:100) {
  
  f <- rnorm(n*k, mean = 0, sd = 1)
  ff <- matrix(f,n,k)
  Y <-  rnorm(n, mean = 0, sd = 1)
  V = as.numeric(Y>0)
  stru1 <- matrix(alpha, n, p-1, byrow = TRUE) + ff %*% Beta + Y %*% t(phi)
  q <- cbind(exp(stru1),1)/(rowSums(exp(stru1))+1)
  X <- matrix(0,n,p)
  for (j in 1:n) {
    count <- sample(1:p, round(runif(1,100,1000)), replace = T, prob = q[j, ])
    for (i in 1:p) {
      X[j,i] <- sum(count==i)
    }
  }
  X.train <- X[1:(n/2),];V.train <- V[1:(n/2)];ff.train <- ff[1:(n/2),]
  X.test <- X[(n/2+1):n,];V.test <- V[(n/2+1):n];ff.test <- ff[(n/2+1):n,]
  Y.train <- Y[1:(n/2)]
  Y.test <- Y[(n/2+1):n]
  
  fit1 <- FAIR(X.train, V.train, penalty = NULL, n.factors=k, maxit = 500)
  phi1 <- fit1[["model.coefs"]][["phi"]]
  Beta1 <- fit1[["model.coefs"]][["Beta"]]
  phi.fit1 <- cbind(phi.fit1, c(phi1,0))
  Beta.fit1 <- cbind(Beta.fit1, c(Beta1,0))
  
  m0 <- rowSums(X.train)
  m1 <- rowSums(X.test)
  z1_test <- tcrossprod(X.test[,-p],phi1)
  z1_test <- as.matrix(z1_test)/m1
  z2_test <- tcrossprod(X.test[,-p],Beta1)
  z2_test <- as.matrix(z2_test)/m1
  Z1 <- cbind(z1_test,z2_test)
  
  z1_train <- tcrossprod(X.train[,-p],phi1)
  z1_train <- as.matrix(z1_train)/m0
  z2_train <- tcrossprod(X.train[,-p],Beta1)
  z2_train <- as.matrix(z2_train)/m0

  
  fit2 <- mnlm(cl=NULL, V.train, X.train,free=1)
  B <- coef(fit2)
  # fit3 <- mnlm(cl, V.train, X.train, cv=T, nfold=10)
  # B <- coef(fit3, select="1se")
  phi2 <- as.vector(B[-1,,drop=FALSE])
  phi.fit2 <- cbind(phi.fit2, phi2)
  z2_test <- X.test%*%phi2/m1
  Z2 <- z2_test
  z2_train <- X.train%*%phi2/m0
  
  ##############
  V.train.new1 <- cbind(ff.train,V.train)
  goldfit <- MGLMreg(formula=X.train~ff.train+V.train, data = data.frame(V.train.new1), dist="MN")
  coef1 <- coef(goldfit)[-1,]
  phi3 <- coef1[2,]
  Beta3 <- coef1[1,]
  phi.fit3 <- cbind(phi.fit3, c(phi3,0))
  Beta.fit3 <- cbind(Beta.fit3, c(Beta3,0))

  z3_test <- X.test[,-p]%*%t(coef1)/m1
  Z3 <- z3_test
  z3_train <- X.train[,-p]%*%t(coef1)/m0

  
  ######################
  spantrue <- cbind(phi,as.vector(Beta))
  spanestim1 <- cbind(as.vector(phi1),as.vector(Beta1))
  spanestim2 <- cbind(phi2-phi2[p],0)[-p,]
  spanestim3 <- t(coef1)
  
  eval1 <- cbind(eval1,eval.space(spanestim1,spantrue))
  eval2 <- cbind(eval2,eval.space(spanestim2,spantrue))
  eval3 <- cbind(eval3,eval.space(spanestim3,spantrue))

  
  
  
  D1 <- c(D1,KCI(X.test/m1, Y.test, Z1)$pvalue)
  D2 <- c(D2,KCI(X.test/m1, Y.test, Z2)$pvalue)
  D3 <- c(D3,KCI(X.test/m1, Y.test, Z3)$pvalue)

  
})

eval1 <- eval1[,-1]
eval2 <- eval2[,-1]
eval3 <- eval3[,-1]

table1data <- data.frame(data=rep(paste("Ex.1 n =",n,"p =",p),3),method=c("FAIR","MNIR","MNIR-f"),
                         VCC.mean=c(mean(eval1[1,]),mean(eval2[1,]),mean(eval3[1,])),
                         VCC.sd=c(sd(eval1[1,]),sd(eval2[1,]),sd(eval3[1,])),
                         TCC.mean=c(mean(eval1[2,]),mean(eval2[2,]),mean(eval3[2,])),
                         TCC.sd=c(sd(eval1[2,]),sd(eval2[1,]),sd(eval3[2,])),
                         D.mean=c(mean(eval1[2,]),mean(eval2[2,]),mean(eval3[2,])),
                         D.sd=c(sd(eval1[2,]),sd(eval2[1,]),sd(eval3[2,])))


table1data
