##k=0
library(MGLM)
source("FAIR.R")

n=60; p=10; k=1;
nfbic <- numeric();nfaic <- numeric();nfhic1 <- numeric();nfhic2 <- numeric()
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
  bic <- numeric();aic <- numeric()
  hic1 <- numeric();hic2 <- numeric()
  for (i in 0:4) {
    res1 <- FAIR(X.train, V.train, base = p, n.factors=i, maxit = 500)
    bic[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn)+(p-1)*k2)*log(nn)
    aic[i+1] <- -2*res1[["ELBO"]] + (i*((p-1)+2*nn)+(p-1)*k2)*2
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
      pd <- 2*logLik(sweep4)-2*(res1[["ELBO1"]]+ sum(lgamma(m.train + 1)) - sum(lgamma(X.train + 1)))
      hic1[i+1] <- -2*logLik(sweep4)+2*(pd+i*(p-1)+(p-1)*k2)
      hic2[i+1] <- -2*logLik(sweep4)+log(nn)*(i*(p-1)+(p-1)*k2+pd)
    }
  }
  nfbic[runs] <- which.min(bic)-1
  nfaic[runs] <- which.min(aic)-1
  nfhic1[runs] <- which.min(hic1)-1
  nfhic2[runs] <- which.min(hic2)-1
}

table2data <- data.frame(data=rep(paste("Ex.1 n =",n,"p =",p),4),method=c("AIC","BIC","HIC_1","HIC_2"),
                         d0=c(sum(nfaic==0),sum(nfbic==0),sum(nfhic1==0),sum(nfhic2==0)),
                         d1=c(sum(nfaic==1),sum(nfbic==1),sum(nfhic1==1),sum(nfhic2==1)),
                         d2=c(sum(nfaic==2),sum(nfbic==2),sum(nfhic1==2),sum(nfhic2==2)),
                         dl2=c(sum(nfaic>2),sum(nfbic>2),sum(nfhic1>2),sum(nfhic2>2)))

table2data



