##k=0
library(MASS)
library(MGLM)
source("FAIR.R")

n=60; p=10; k=1;
nf11 <- numeric();nf12 <- numeric();nf13 <- numeric();nf14 <- numeric();nf15 <- numeric();
nf16 <- numeric();nf17 <- numeric();nf18 <- numeric();nf19 <- numeric();nf20 <- numeric();
set.seed(927)
alpha <- rnorm(p-1, mean = 0, sd = 1)
# phi <- sample(c(-3:-1,1:3),p-1,replace = T)
phi <- sample(c(-1:1),p-1,replace = T)
bb.seq <- -1:1#c(seq(from=-3,to=-1,by=1),seq(from=1,to=3,by=1))
Beta <- matrix(sample(bb.seq,k*(p-1),replace = TRUE),k,p-1)
alpha <- rep(0,p-1)
phi <- c(1,1,1,1,0,0,0,0,0)

# Beta <- t(cbind(c(-1,1,-1,1,0,0,0,0,0),c(0,0,0,0,1,1,1,1,1)))
Beta <- c(-1,1,-1,1,1,1,1,1,1)

for (runs in 1:100) {
  f <- rnorm(n*k, mean = 0, sd = 1)
  ff <- matrix(f,n,k)
  Y <- rnorm(n, mean = 0, sd = 1)
  V = Y
  stru1 <- matrix(alpha, n, p-1, byrow = TRUE) + ff %*% Beta + V %*% t(phi)
  q <- cbind(exp(stru1),1)/(rowSums(exp(stru1))+1)
  
  X <- matrix(0,n,p)
  for (j in 1:n) {
    count <- sample(1:p, round(runif(1,1000,1e5)), replace = T, prob = q[j, ])
    for (i in 1:p) {
      X[j,i] <- sum(count==i)
    }
  }
  
  X.train <- X[1:(n/2),];V.train <- V[1:(n/2)];ff.train <- ff[1:(n/2),]
  X.test <- X[(n/2+1):n,];V.test <- V[(n/2+1):n];ff.test <- ff[(n/2+1):n,]
  m.train=rowSums(X.train)
  k2 <- ncol(as.matrix(V.train))
  nn <- nrow(as.matrix(V.train))
  bic10 <- numeric();bic11 <- numeric();bic12 <- numeric();aic11 <- numeric();aic12 <- numeric();
  bic20 <- numeric();bic21 <- numeric();bic22 <- numeric();aic21 <- numeric();aic22 <- numeric();
  hic <- numeric();dic1 <- numeric();dic2 <- numeric();dic3 <- numeric();dic4 <- numeric();
  for (i in 0:4) {
    res1 <- FAIR(X.train, V.train, base = p, n.factors=i, maxit = 500)
    bic10[i+1] <- -2*res1[["ELBO"]]
    bic11[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn))*log(nn)
    bic12[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn)+(p-1)*k2)*log(nn)
    aic11[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn))*2
    aic12[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn)+(p-1)*k2)*2
    if(i==0){
      V.new=as.matrix(V.train)
    }else{
      new.va.mu <- res1[["va.coefs"]][["Mu"]]
      f.hat <- matrix(new.va.mu, nn, i)
      # ff.hat <- t(t(f.hat)/sqrt(colSums(f.hat^2)))*sqrt(nn)
      V.new <- cbind(V.train,f.hat)
    }
    qt=try(sweep4 <- MGLMreg(formula=X.train~1+V.new, dist="MN"))
    if("try-error" %in% class(qt)){ print("error")
    }else{
      hic[i+1] <- -2*logLik(sweep4)+log(nn)*(i*(p-1)+(p-1)*k2)+2*i*nn
      pd <- 2*logLik(sweep4)-2*(res1[["ELBO1"]]+ sum(lgamma(m.train + 1)) - sum(lgamma(X.train + 1)))
      dic1[i+1] <- -2*logLik(sweep4)+2*(pd+i*(p-1)+(p-1)*k2)
      dic2[i+1] <- -2*logLik(sweep4)+log(nn)*(i*(p-1)+(p-1)*k2)+2*pd
      dic3[i+1] <- -2*logLik(sweep4)+2*(i*(p-1)+(p-1)*k2)+log(nn)*pd
      dic4[i+1] <- -2*logLik(sweep4)+log(nn)*(i*(p-1)+(p-1)*k2+pd)
    }
  }
  nf11[runs] <- which.min(bic10)-1
  nf12[runs] <- which.min(bic11)-1
  nf13[runs] <- which.min(bic12)-1
  nf14[runs] <- which.min(aic11)-1
  nf15[runs] <- which.min(aic12)-1
  nf16[runs] <- which.min(hic)-1
  nf17[runs] <- which.min(dic1)-1
  nf18[runs] <- which.min(dic2)-1
  nf19[runs] <- which.min(dic3)-1
  nf20[runs] <- which.min(dic4)-1
}

sum(nf11==0)
sum(nf12==0)
sum(nf13==0)
sum(nf14==0)
sum(nf15==0)
sum(nf16==0)
sum(nf17==0)
sum(nf18==0)
sum(nf19==0)
sum(nf20==0)

sum(nf11==1)
sum(nf12==1)
sum(nf13==1)
sum(nf14==1)
sum(nf15==1)
sum(nf16==1)
sum(nf17==1)
sum(nf18==1)
sum(nf19==1)
sum(nf20==1)


sum(nf11==2)
sum(nf12==2)
sum(nf13==2)
sum(nf14==2)
sum(nf15==2)
sum(nf16==2)
sum(nf17==2)
sum(nf18==2)
sum(nf19==2)
sum(nf20==2)

sum(nf11>2)
sum(nf12>2)
sum(nf13>2)
sum(nf14>2)
sum(nf15>2)
sum(nf16>2)
sum(nf17>2)
sum(nf18>2)
sum(nf19>2)
sum(nf20>2)

save.image("nfEx1.RData")



