# Code for Factor augmented inverse regression
# Author: Daolin Pang
#============================================================##

FAIR <- function(X, Y, base = NULL, penalty=NULL, trace = FALSE, n.factors=2,ngridpt=60,keep.path=FALSE, maxit = 1000) {
  n <- nrow(X)
  p <- ncol(X)-1
  k <- n.factors
  if (is.null(base)) base=p+1
  if (base != (p+1)&base!=1){
    X <- X[,c(1:(base-1),(base+1):(p+1),base)]
  }
  if(base==1){
    X <- X[,c((base+1):(p+1),base)]
  }
  m <- rowSums(X)
  V <- as.matrix(Y)
  k2 <- ncol(V)
  ###Initialization###
  
  # alpha <- rep(1,p)
  # phi <- rep(1,p)
  # Beta <- rep(1, k*p)
  # alpha0 <- c(0,1,2,3,1.5)
  # phi0 <- c(-1,0,1,2,2)
  # Beta0 <- matrix(rep(c(1,2,3,0.5,2), 2), 2)
  alpha0 <- rep(1, p)
  phi0 <- rep(1, k2*p)
  Beta0 <- rep(1, k*p)
  new.model.coefs <- c(phi0, alpha0, Beta0)
  new.va.mu <- rep(1, n*k)
  new.va.sigma <- rep(1, n*k)
  
  ###obj fun and grad fun###
  obj <- function(model.coefs, va.mu = NULL, va.sigma = NULL){
    
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru1 <- matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi
    stru2 <- exp(stru1 + 0.5 * Sigma %*% Beta^2)
    
    y <- sum(X[,1:p] * stru1) - sum(m * log(rowSums(stru2)+1)) - 0.5 * sum(Mu^2 + Sigma) + 0.5 * sum(log(Sigma)) + n*k/2
    
    return(y)
  }
  
  obj_model_coefs <- function(model.coefs, va.mu = NULL, va.sigma = NULL){
    
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru1 <- matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi
    stru2 <- exp(stru1 + 0.5 * Sigma %*% Beta^2)
    
    y <- sum(X[,1:p] * stru1) - sum(m*log(rowSums(stru2)+1))
    
    return(y)
  }
  
  gr_model_coefs <- function(model.coefs, va.mu = NULL, va.sigma = NULL) {
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru1 <- exp(matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi + 0.5 * Sigma %*% Beta^2)
    stru2 <- m/(rowSums(stru1)+1)*stru1
    
    ###gr_alpha
    y1 <- colSums(X[,1:p]) - colSums(stru2)
    ###gr_phi
    y2 <- crossprod(V, X[,1:p]) - crossprod(V, stru2)
    ###gr_Beta
    y3 <- crossprod(Mu, X[,1:p]) - crossprod(Mu, stru2) - crossprod(Sigma, stru2) * Beta
    
    
    yy <- c(c(y1), c(y2), c(y3))
    return(yy)
  }
  
  obj_Mu <- function(model.coefs, va.mu = NULL, va.sigma = NULL){
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru2 <- exp(matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi + 0.5 * Sigma %*% Beta^2)
    
    y <- sum(X[,1:p] * Mu %*% Beta) - sum(m * log(rowSums(stru2)+1)) - 0.5 * sum(Mu^2)
    
    return(y)
  }
  
  gr_Mu <- function(model.coefs, va.mu = NULL, va.sigma = NULL) {
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru1 <- exp(matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi + 0.5 * Sigma %*% Beta^2)
    y <- tcrossprod(X[,1:p], Beta) - m/(rowSums(stru1)+1) * tcrossprod(stru1, Beta) - Mu
    
    return(c(y))
  }
  
  obj_Sigma <- function(model.coefs, va.mu = NULL, va.sigma = NULL){
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru2 <- exp(matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi + 0.5 * Sigma %*% Beta^2)
    
    y <- -sum(m * log(rowSums(stru2)+1)) - 0.5 * sum(Sigma) + 0.5 * sum(log(Sigma))
    
    return(y)
  }
  
  gr_Sigma <- function(model.coefs, va.mu = NULL, va.sigma = NULL) {
    alpha <- model.coefs[1:p]
    phi <- matrix(model.coefs[(p+1):(p+k2*p)], k2, p)
    Beta <- matrix(model.coefs[-(1:((k2+1)*p))], k, p)
    Mu <- matrix(va.mu, n, k)
    Sigma <- matrix(va.sigma^2, n, k)
    
    stru1 <- exp(matrix(alpha, n, p, byrow = TRUE) + Mu %*% Beta + V %*% phi + 0.5 * Sigma %*% Beta^2)
    
    # y <- matrix(abs(va.sigma), n, k)*(-m/rowSums(stru1)*tcrossprod(stru1, Beta^2) - 1 + 1 / Sigma)
    y <- matrix(va.sigma, n, k)*(-m/(rowSums(stru1)+1)*tcrossprod(stru1, Beta^2) - 1 + 1 / Sigma)
    
    return(c(y))
  }
  
  ###VA iteration
  cur.ELBO <- -1e9; iter <- 1; ratio <- 10; diff=1e5; eps = 1e-5; max.iter = maxit;
  mu.cur.logfunc <- -1e9;sigma.cur.logfunc <- -1e9;model.coefs.cur.logfunc <- -1e9;
  
  while((diff> eps*(abs(cur.ELBO)+eps)) && iter <= max.iter) {
    if(trace) cat("Iteration:", iter, "\n")
    
    ###optim mu
    q <- try(optim(new.va.mu, model.coefs=new.model.coefs, va.sigma=new.va.sigma, method = "BFGS", fn = obj_Mu, gr = gr_Mu, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && mu.cur.logfunc > q$value){if(trace)
        cat("Optimization of mu did not improve on iteration step ",iter,"\n")
      }else{
        if(trace) cat("Variational parameters mu updated","\n")
        new.va.mu <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of mu did not converge on iteration step ", iter,"\n") }
      }
    }
    sigma.cur.logfunc <- obj_Sigma(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ###optim sigma
    q <- try(optim(new.va.sigma, model.coefs=new.model.coefs, va.mu=new.va.mu, method = "BFGS", fn = obj_Sigma, gr = gr_Sigma, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && sigma.cur.logfunc > q$value){if(trace)
        cat("Optimization of sigma did not improve on iteration step ",iter,"\n")
      }else{
        if(trace) cat("Variational parameters sigma updated","\n")
        new.va.sigma <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of sigma did not converge on iteration step ", iter,"\n") }
      }
    }
    model.coefs.cur.logfunc <- obj_model_coefs(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ###optim model.coefs
    q <- try(optim(new.model.coefs, va.mu=new.va.mu, va.sigma=new.va.sigma, method = "BFGS", fn = obj_model_coefs, gr = gr_model_coefs, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && model.coefs.cur.logfunc > q$value){if(trace)
        cat("Optimization of model.coefs did not improve on iteration step ",iter,"\n");
      }else{
        if(trace) cat("Variational parameters model.coefs updated","\n")
        new.model.coefs <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of model.coefs did not converge on iteration step ", iter,"\n") }
      }
    }
    
    
    mu.cur.logfunc <- obj_Mu(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ## Take values of ELBO to define stopping rule
    new.ELBO <- obj(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    diff=abs(new.ELBO-cur.ELBO)
    ratio <- abs(new.ELBO/cur.ELBO);
    if(trace) cat("New ELBO:", new.ELBO,"cur ELBO:", cur.ELBO, "Ratio of ELBO", ratio, ". Difference in ELBO:",diff,"\n")
    cur.ELBO <- new.ELBO
    iter <- iter + 1
    
  }
  
  new.model.coefs[-(1:((k2+1)*p))] <- Beta0
  new.va.mu <- rep(1, n*k)
  new.va.sigma <- rep(1, n*k)
  
  cur.ELBO <- -1e9; iter <- 1; ratio <- 10; diff=1e5; eps = 1e-5; max.iter = maxit;
  mu.cur.logfunc <- -1e9;sigma.cur.logfunc <- -1e9;model.coefs.cur.logfunc <- -1e9;
  
  while((diff> eps*(abs(cur.ELBO)+eps)) && iter <= max.iter) {
    if(trace) cat("Iteration:", iter, "\n")
    
    ###optim mu
    q <- try(optim(new.va.mu, model.coefs=new.model.coefs, va.sigma=new.va.sigma, method = "BFGS", fn = obj_Mu, gr = gr_Mu, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && mu.cur.logfunc > q$value){if(trace)
        cat("Optimization of mu did not improve on iteration step ",iter,"\n")
      }else{
        if(trace) cat("Variational parameters mu updated","\n")
        new.va.mu <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of mu did not converge on iteration step ", iter,"\n") }
      }
    }
    sigma.cur.logfunc <- obj_Sigma(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ###optim sigma
    q <- try(optim(new.va.sigma, model.coefs=new.model.coefs, va.mu=new.va.mu, method = "BFGS", fn = obj_Sigma, gr = gr_Sigma, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && sigma.cur.logfunc > q$value){if(trace)
        cat("Optimization of sigma did not improve on iteration step ",iter,"\n")
      }else{
        if(trace) cat("Variational parameters sigma updated","\n")
        new.va.sigma <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of sigma did not converge on iteration step ", iter,"\n") }
      }
    }
    model.coefs.cur.logfunc <- obj_model_coefs(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ###optim model.coefs
    q <- try(optim(new.model.coefs, va.mu=new.va.mu, va.sigma=new.va.sigma, method = "BFGS", fn = obj_model_coefs, gr = gr_model_coefs, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
    if("try-error" %in% class(q)){ print("error")
    }else{
      if(iter > 1 && model.coefs.cur.logfunc > q$value){if(trace)
        cat("Optimization of model.coefs did not improve on iteration step ",iter,"\n");
      }else{
        if(trace) cat("Variational parameters model.coefs updated","\n")
        new.model.coefs <- q$par;
        if(q$convergence != 0) { if(trace) cat("Optimization of model.coefs did not converge on iteration step ", iter,"\n") }
      }
    }
    
    
    
    mu.cur.logfunc <- obj_Mu(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    
    ## Take values of ELBO to define stopping rule
    new.ELBO <- obj(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
    diff=abs(new.ELBO-cur.ELBO)
    ratio <- abs(new.ELBO/cur.ELBO);
    if(trace) cat("New ELBO:", new.ELBO,"cur ELBO:", cur.ELBO, "Ratio of ELBO", ratio, ". Difference in ELBO:",diff,"\n")
    cur.ELBO <- new.ELBO
    iter <- iter + 1
    
  }
  
  if(iter > max.iter){
    print("LMNIR Not converging!")
  }
  f.hat <- matrix(new.va.mu, n, k)
  ff.hat <- t(t(f.hat)/sqrt(colSums(f.hat^2)))*sqrt(n)
  
  out.list <- list()
  out.list$ELBO <- cur.ELBO
  out.list$ELBO1 <- obj_model_coefs(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = new.va.sigma)
  out.list$lc1 <- obj_model_coefs(model.coefs = new.model.coefs, va.mu = new.va.mu, va.sigma = 0)
  
  out.list$n.factor <- n.factors
  out.list$iter=iter-1
  out.list$model.coefs$alpha <- new.model.coefs[1:p]
  out.list$model.coefs$phi <- matrix(new.model.coefs[(p+1):(p+k2*p)], k2, p)
  out.list$model.coefs$Beta <- matrix(new.model.coefs[-(1:((k2+1)*p))], k, p) 
  
  if(!is.null(penalty)){
    
    V.new <- cbind(V,f.hat)
    V.new <- t(t(V.new)/sqrt(colSums(V.new^2)))*sqrt(n)
    
    sweep4 <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty="group", va.sigma=0, ngridpt=ngridpt,keep.path = keep.path, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    
    new.model.coefs4 <- c(t(sweep4@select@coefficients))
    out.list$model.coefs4$alpha <- new.model.coefs4[1:p]
    out.list$model.coefs4$phi <- matrix(new.model.coefs4[(p+1):(p+k2*p)], k2, p,byrow = T)
    out.list$model.coefs4$Beta <- matrix(new.model.coefs4[-(1:((k2+1)*p))], k, p,byrow = T)
    out.list$groupsolution <- sweep4
    new.model.coefs4[-(1:((k2+1)*p))] <- as.vector(matrix(new.model.coefs4[-(1:((k2+1)*p))], k, p,byrow = T))
    new.va.mu4 <- t(t(f.hat)/sqrt(colSums(f.hat^2)))*sqrt(n)
    new.va.sigma4 <- matrix(new.va.sigma, n, k)
    new.va.sigma4 <- t(t(new.va.sigma4)/sqrt(colSums(f.hat^2)))*sqrt(n)
    V <- t(t(V)/sqrt(colSums(V^2)))*sqrt(n)
    out.list$ELBO2 <- obj_model_coefs(model.coefs = new.model.coefs4, va.mu = as.vector(new.va.mu4), va.sigma = as.vector(new.va.sigma4))
    out.list$lc2 <- obj_model_coefs(model.coefs = new.model.coefs4, va.mu = as.vector(new.va.mu4), va.sigma = 0)
    
  }
  
  ## print the output
  
  out.list$va.coefs$Mu <- ff.hat
  # out.list$va.coefs$Sigma <- matrix(new.va.sigma^2, n, k)
  
  
  return(out.list)
}

LFIRtune <- function(formula, data, k, dist, va.sigma, penalty, lambdas, tau, ngridpt, warm.start = TRUE, 
                     keep.path = FALSE, display = FALSE, init, penidx, ridgedelta, maxiters = 150, 
                     epsilon = 1e-05, regBeta = FALSE, overdisp) {
  
  call <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame(n = 1))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  X <- model.matrix(mt, mf, contrasts)
  ow <- getOption("warn")
  
  if (missing(tau)) {
    if (penalty == "SCAD") {
      tau=3.7
    } else {
      tau=3
    }
  }
  # if (!missing(weight)) 
  #   weight <- weight[rowSums(Y) != 0]
  X <- as.matrix(X[rowSums(Y) != 0, ])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
  
  d <- ncol(Y)
  m <- rowSums(Y)
  p <- ncol(X)
  N <- nrow(Y)
  # if (missing(weight)) 
  #   weight <- rep(1, N)
  ## ----------------------------------------## 
  ## Check distribution and d
  ## ----------------------------------------## 
  if (dist == "GDM" && d == 2) 
    stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
  ## ----------------------------------------## 
  ## Pennalty type
  ## ----------------------------------------## 
  if (penalty == "group")
    penalty <- "group_row"
  # if (!penalty %in% c("sweep", "group_row", "nuclear")) 
  #   stop("Penalty type can only be sweep, group, or nuclear.")
  ## ----------------------------------------##
  ## regularization set
  ## ----------------------------------------##
  if (missing(penidx)) 
    penidx <- rep(TRUE, p)
  ## ----------------------------------------## 
  ## Starting point
  ## ----------------------------------------## 
  if (missing(init)) {
    if (dist == "MN") 
      init <- matrix(0, p, (d - 1)) else if (dist == "DM") 
        init <- matrix(0, p, d) else if (dist == "GDM") 
          init <- matrix(0, p, 2 * (d - 1)) else if (dist == "NegMN") {
            if (regBeta) 
              init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
          }
  }
  if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
    options(warn = -1)
    est <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, init = init, X = X, weight = weight, 
                     epsilon = epsilon, maxiters = maxiters, parallel = FALSE, cores = 1, 
                     cl = NULL, sys = NULL, display = FALSE))
    overdisp <- est$coefficients$phi
    options(warn = ow)
  } else {
    overdisp <- NULL
  }
  ## ----------------------------------------## 
  ## Ridgedelta
  ## ----------------------------------------## 
  if (missing(ridgedelta)) 
    ridgedelta <- 1/(p * d)
  ## ----------------------------------------## 
  ## Pennalty values
  ## ----------------------------------------## 
  fit.max <- NULL
  if (missing(lambdas)) {
    if (missing(ngridpt)) 
      ngridpt <- 10
    ## ----------------------------------------## 
    ## Find maximum lambda
    ## ----------------------------------------## 
    fit.max <- eval(call("LFIRlasso.fit", Y = Y, X = X, k = k, dist = dist, va.sigma = va.sigma, penalty = penalty, lambda = Inf, 
                         tau = tau, penidx = penidx, init = init, ridgedelta = ridgedelta, 
                         maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
    
    maxlambda <- fit.max@maxlambda
    # lambdas <- exp(seq(from=log(maxlambda/N), to=log(maxlambda), length.out=15))
    # ----------------------------------------## Find minimum lambda
    # if(penalty=='group_row'){ for(j in 1:15){ if(j==1) B0 <- fit.max$coefficients
    # else B0 <- B_hat temp <- eval(call('MGLMsparsereg.fit', Y=Y, X=X, dist=dist,
    # lambda=lambdas[j], penalty=penalty, weight=weight, penidx=penidx, init=init,
    # ridgedelta=ridgedelta, maxiters=maxiters, epsilon=epsilon, regBeta=regBeta,
    # overdisp=overdisp)) B_hat <- temp$coefficients nz <- sum(rowSums(B_hat^2)>0)
    # if(nz==p){ next }else{ if(j>1) minlambda <- lambdas[j-1] else minlambda <-
    # lambdas[1]/10 break } } }else{
    minlambda <- maxlambda/N
    # minlambda <- maxlambda/N/2#0617
    # minlambda <- 1
    # }
    lambdas <- exp(seq(from = log(maxlambda), to = log(minlambda), length.out = ngridpt))
  } else {
    ngridpt <- length(lambdas)
  }
  
  BICs <- rep(NA, ngridpt)
  AICs <- rep(NA, ngridpt)
  logL <- rep(NA, ngridpt)
  Dof <- rep(NA, ngridpt)
  select.list <- list()
  
  for (j in 1:ngridpt) {
    if (j == 1 & !is.null(fit.max)) {
      temp <- fit.max
      select.list[[j]] <- temp
      B_hat <- temp@coefficients
      BICs[j] <- temp@BIC
      AICs[j] <- temp@AIC
      logL[j] <- temp@logL
      Dof[j] <- temp@Dof
      if (display) 
        print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", 
                    sprintf("%.2f", BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", 
                    sprintf("%.2f", logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
      next
    }
    
    if (warm.start) {
      if (j == 1) 
        B0 <- init else B0 <- B_hat
    } else {
      B0 <- init
    }
    temp <- eval(call("LFIRlasso.fit", Y = Y, X = X, k = k, dist = dist, va.sigma, penalty = penalty, lambda = lambdas[j], 
                      tau = tau, penidx = penidx, init = B0, ridgedelta = ridgedelta, 
                      maxiters = maxiters, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
    select.list[[j]] <- temp
    B_hat <- temp@coefficients
    BICs[j] <- temp@BIC
    AICs[j] <- temp@AIC
    logL[j] <- temp@logL
    Dof[j] <- temp@Dof
    
    if (display) 
      print(paste(j, " lamda=", sprintf("%.2f", lambdas[j]), "  BIC=", sprintf("%.2f", 
                                                                               BICs[j]), " AIC=", sprintf("%.2f", AICs[j]), " logL=", sprintf("%.2f", 
                                                                                                                                              logL[j]), " Dof=", sprintf("%.2f", Dof[j]), sep = ""))
    
  }
  chosen.lambda <- lambdas[which.min(BICs)]
  select <- select.list[[which.min(BICs)]]
  select@call <- match.call()
  select@data <- list(Y = Y, X = X)
  select@distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
                                                                    "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
                                                                                                    "Negative Multinomial")))
  select@lambda <- chosen.lambda
  
  outDf <- data.frame(Dof = Dof, Lambda = lambdas, BIC = BICs, AIC = AICs, logL = logL)
  
  
  ## ----------------------------------------## 
  ## Set class
  ## ----------------------------------------## 
  if (keep.path) new("LFIRtune", call = select@call, select = select, path = outDf, select.list = select.list)
  else new("LFIRtune", call = select@call, select = select, path = outDf, select.list = list())
  
} 

LFIRlasso <- function(formula, data, k, dist, va.sigma, penalty, lambda, tau, init, penidx, 
                      maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
  call <- match.call()
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  # m <- match(c("formula", "data", "weight"), names(mf), 0L)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame(n = 1))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  X <- model.matrix(mt, mf, contrasts)
  # if (!missing(weight)) 
  #   weight <- weight[rowSums(Y) != 0]
  X <- as.matrix(X[rowSums(Y) != 0, ])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
  d <- ncol(Y)
  m <- rowSums(Y)
  p <- ncol(X)
  N <- nrow(Y)
  if (penalty == "group")
    penalty <- "group_row"
  if (missing(tau)) {
    if (penalty == "SCAD") {
      tau=3.7
    } else {
      tau=3
    }
  }
  
  if (missing(penidx)) 
    penidx <- rep(TRUE, p)
  if (missing(init)) {
    if (dist == "MN") {
      init <- matrix(0, p, (d - 1))
    } else if (dist == "DM") {
      init <- matrix(0, p, d)
    } else if (dist == "GDM") {
      init <- matrix(0, p, 2 * (d - 1))
    } else if (dist == "NegMN") {
      if (regBeta) 
        init <- matrix(0, p, (d + 1)) else init <- matrix(0, p, d)
    }
  }
  if (dist == "NegMN" && regBeta == FALSE && missing(overdisp)) {
    est <- DMD.NegMN.fit(Y)
    overdisp <- est$estimate[d + 1]
  } else {
    overdisp <- NULL
  }
  if (missing(ridgedelta)) 
    ridgedelta <- 1/(p * d)
  est <- eval(call("LFIRlasso.fit", Y = Y, X = X, k = k, dist = dist, va.sigma = va.sigma, penalty = penalty,  
                   lambda = lambda, tau = tau, init = init, penidx = penidx, maxiters = maxiters, 
                   ridgedelta = ridgedelta, epsilon = epsilon, regBeta = regBeta, overdisp = overdisp))
  est@call <- match.call()
  
  ## ----------------------------------------## 
  ## Set class
  ## ----------------------------------------## 
  return(est) 
  
}

## ============================================================## 
## The FAIR sparse reg function 
## ============================================================##
LFIRlasso.fit <- function(Y, X, k, dist, va.sigma, penalty, lambda, tau, init, penidx, 
                          maxiters = 150, ridgedelta, epsilon = 1e-05, regBeta = FALSE, overdisp) {
  
  tol <- 1e-8 
  d <- ncol(Y)
  m <- rowSums(Y)
  p <- ncol(X)
  N <- nrow(Y)
  
  
  beta_old <- init
  B <- init
  alpha <- 1
  alpha_iter <- list()
  objval <- Inf
  niter <- 1
  isdescent <- TRUE
  oriRidge <- ridgedelta
  while ((niter < 3) || (!stop)) {
    niter <- niter + 1
    beta_old <- B
    obj1 <- objval
    if (niter <= 2) 
      S <- init
    loss <- LFIR.loss(Y, X, k, S, dist, regBeta, overdisp)
    loss.S <- loss[[1]]
    loss.D1S <- loss[[2]]
    for (l in 1:50) {
      A <- S + ridgedelta * loss.D1S#- ridgedelta * loss.D1S #0617
      B <- matrix(NA, nrow(B), ncol(B))
      B[!penidx, ] <- A[!penidx, ]
      if (d > 2) {
        Apen <- A[penidx, ]
      } else if (d == 2) 
        Apen <- matrix(A[penidx, ], , 1)
      pen <- matrix_threshold(X = Apen, penalty = penalty, lambda = ridgedelta * lambda, tau = tau)
      B[penidx, ] <- pen[[1]]
      penval <- pen[[2]]
      if (all(abs(B[penidx, ]) < tol)) {
        if (penalty == "group_row") {
          maxlambda <- max(sqrt(rowSums(Apen^2)))/ridgedelta
        } else {
          maxlambda <- max(abs(Apen))/ridgedelta
        }
        # if (penalty == "sweep") {
        # } else  if (penalty == "nuclear") {
        #   sing.vec <- svd(Apen)$d
        #   maxlambda <- sing.vec[1]/ridgedelta
        # }
        if (lambda > maxlambda) {
          B[penidx, ] <- 0
          penval <- 0
        }
      } else {
        maxlambda <- numeric() # 12/1/17
      }
      objloss <- LFIR.loss(Y, X, k, B, dist, va.sigma, regBeta, overdisp)
      loss.S2 <- objloss[[1]]
      loss.D1S2 <- objloss[[2]]
      objval <- loss.S2 + penval
      BminusS <- B - S
      surval <- loss.S - sum(loss.D1S * BminusS) + norm(BminusS, type = "F")^2/2/ridgedelta + 
        penval#+ sum(loss.D1S * BminusS)
      if (!is.na(objval) && !is.na(surval)) 
        if (objval <= surval) 
          break else ridgedelta <- ridgedelta/2
    }
    alpha_old <- alpha
    alpha <- (1 + sqrt(4 + alpha_old^2))/2
    if (!is.na(objval) & objval <= obj1) {
      stop <- abs(obj1 - objval) < epsilon * (abs(obj1) + 1)
      S <- B + (alpha_old - 1)/alpha * (B - beta_old)
    } else {
      objval <- obj1
      if (isdescent) {
        isdescent <- FALSE
        stop <- FALSE
        S <- B + (alpha_old - 1)/alpha * (beta_old - B)
        B <- beta_old
      } else {
        stop <- TRUE
      }
    }
    if (niter >= maxiters) 
      stop <- TRUE
  }
  if (d > 2) 
    Bpen <- B[penidx, ] else if (d == 2) 
      Bpen <- matrix(B[penidx, ], , 1)
  if (penalty == "group_row") {
    Dof <- sum(!penidx) * ncol(B) + sum(colSums(Bpen^2) != 0) + sum((nrow(B) -
                                                                       1) * colSums(Bpen^2)/colSums(Apen^2))
  } else {
    Dof <- sum(!penidx) * ncol(B) + sum(Bpen != 0)
  }
  
  logL <- -LFIR.loss(Y, X, k, B, dist, va.sigma, regBeta, overdisp)[[1]]
  AIC <- -2 * logL + 2 * Dof
  BIC <- -2 * logL + log(N) * Dof
  
  distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
                                                             "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
                                                                                             "Negative Multinomial")))
  # penalty_name <- ifelse(penalty == "group_row", "group", ifelse(penalty == "nuclear", 
  #                                                                "nuclear", "sweep"))
  new("LFIRlasso", call = match.call(), data = list(Y = Y, X = X), 
      coefficients = B, logL = logL, BIC = BIC, AIC = AIC, Dof = Dof, 
      iter = niter, maxlambda = maxlambda, lambda = lambda, 
      distribution = distribution) # , Beta = est$Beta) 
}

## ============================================================================##
## LFIR loss
## ============================================================================##
LFIR.loss <- function(Y, X, k, beta, dist, va.sigma, regBeta = FALSE) {
  
  # if (missing(weight))
  #   weight <- rep(1, nrow(Y))
  
  losslist <- switch(dist,
                     "MN"= ELBO.MN(Y, X, k, beta, dist, va.sigma))
  
  return(losslist)
}


## ============================================================================##
## matrix thresholding
## ============================================================================##

matrix_threshold <- function(X, penalty, lambda, tau) {
  
  N <- nrow(X)
  d <- ncol(X)
  if (penalty == "lasso") {
    B <- lasso_thresholding(X, lambda)
    penalty_value <- lambda * sum(abs(B))
  } else if (penalty == "SCAD") {
    SCAD.pen <- SCAD_thresholding(X, lambda, tau)
    B <- SCAD.pen[[1]]
    penalty_value <- SCAD.pen[[2]]
  } else if (penalty == "MCP") {
    MCP.pen <- MCP_thresholding(X, lambda, tau)
    B <- MCP.pen[[1]]
    penalty_value <- MCP.pen[[2]]
  } else if (penalty == "group_row" || penalty == "group") {
    row_12norm <- sqrt(colSums(X^2))
    vec <- 1 - lambda/row_12norm
    vec[vec < 0] <- 0
    vec <- matrix(vec,nrow = nrow(X),ncol = ncol(X),byrow = T)
    B <- X * vec
    penalty_value <- lambda * sum(sqrt(colSums(B^2)))
  }
  
  return(list(B, penalty_value))
}

## ============================================================================##
## lsq_threshold We can only work on lasso for now
## ============================================================================##
lasso_thresholding <- function(b, lambda) {
  if (lambda < 0) 
    stop("Penalty constant lambda should be nonnegative.")
  
  B <- b
  # print(b)
  B[abs(b) <= lambda] <- 0
  B[b > lambda] <- B[b > lambda] - lambda
  B[b < -lambda] <- B[b < -lambda] + lambda
  
  return(B)
}

SCAD_thresholding <- function(b, lambda, tau) {
  if (lambda < 0) 
    stop("Penalty constant lambda should be nonnegative.")
  if (tau <= 2) 
    stop("Penalty (SCAD) constant tau should be greater than 2.")
  B <- b
  B[abs(b) <= lambda] <- 0
  B[b > lambda & b <= 2*lambda] <- B[b > lambda & b <= 2*lambda] - lambda
  B[b < -lambda & b >= -2*lambda] <- B[b < -lambda & b >= -2*lambda] + lambda
  B[abs(b) >= 2*lambda & abs(b) <= tau*lambda] <- 
    sign(B[abs(b) >= 2*lambda & abs(b) <= tau*lambda])*((tau-1)*abs(B[abs(b) >= 2*lambda & abs(b) <= tau*lambda])-tau*lambda)/(tau-2)
  # B[abs(b) >= tau*lambda] <- B[abs(b) >= tau*lambda]
  ab.vec = abs(B)
  tem0.vec = (lambda*ab.vec) * (ab.vec<lambda)
  tem1.vec = ((tau*lambda*(ab.vec-lambda)-(ab.vec^2-lambda^2)/2)/(tau-1)+lambda^2) * (ab.vec>=lambda)*(ab.vec<tau*lambda)
  tem2.vec = ((tau+1)*lambda^2/2) * (ab.vec>=tau*lambda)
  pen.value <- sum(tem0.vec+tem1.vec+tem2.vec)
  return(list(B,pen.value))
}

MCP_thresholding <- function(b, lambda, tau) {
  if (lambda < 0) 
    stop("Penalty constant lambda should be nonnegative.")
  if (tau <= 1) 
    stop("Penalty (MCP) constant tau should be greater than 1.")
  B <- b
  B[abs(b) <= lambda] <- 0
  B[abs(b) > lambda & abs(b) <= tau*lambda] <- sign(B[abs(b) > lambda & abs(b) <= tau*lambda])*(abs(B[abs(b) > lambda & abs(b) <= tau*lambda])-lambda)/(1-1/tau)
  ab.vec = abs(B)
  tem0.vec = (lambda*ab.vec-ab.vec^2/2/tau) * (ab.vec<tau*lambda)
  tem1.vec = (tau*lambda^2/2) * (ab.vec>=tau*lambda)
  pen.value <- sum(tem0.vec+tem1.vec)
  return(list(B,pen.value))
}
##============================================================## 
## Function: loss
##============================================================##
ELBO.MN <- function(Y, X, k, model.coefs, dist, va.sigma){
  n <- nrow(Y)
  p <- ncol(Y)-1
  m <- rowSums(Y)
  k2 <- nrow(model.coefs)-k-1
  
  # alpha <- model.coefs[1,]
  # phi <- model.coefs[(2:(k1+1)),]
  Beta <- matrix(model.coefs[-(1:(k2+1)),],k,p)
  
  Sigma <- matrix(va.sigma^2, n, k)
  
  
  V <- matrix(X[,(2:(k2+1))],n,k2)
  Mu <- matrix(X[,-(1:(k2+1))],n,k)
  
  stru1 <- X %*% model.coefs
  stru2 <- exp(stru1 + 0.5 * Sigma %*% Beta^2)
  stru3 <- m/(rowSums(stru2)+1)*stru2
  # print(model.coefs)
  ###ELBO
  y <- sum(Y[,1:p] * stru1) - sum(m * log(rowSums(stru2)+1))
  ###gr_alpha
  y1 <- colSums(Y[,1:p]) - colSums(stru3)
  ###gr_phi
  y2 <- crossprod(V, Y[,1:p]) - crossprod(V, stru3)
  
  ###gr_Beta
  y3 <- crossprod(Mu, Y[,1:p]) - crossprod(Mu, stru3) - crossprod(Sigma, stru3) * Beta
  
  y1 <- matrix(y1,1,p)
  y2 <- matrix(y2,k2,p)
  y3 <- matrix(y3,k,p)
  yy <- rbind(y1, y2, y3)
  
  # print(y)
  
  # yy <- yy[,-p]
  return(list(-y, yy))
  
}

MN.loss <- function(Y, X, beta, dist, weight){
  
  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  
  P <- matrix(NA, N, d)
  P[, d] <- rep(1, N)
  P[, 1:(d - 1)] <- exp(X %*% beta)
  P <- P/rowSums(P)
  loss <- -sum(weight * MGLM::dmn(Y, P))
  # print(X)
  
  kr1 <- Y[, -d] - P[, 1:(d - 1)] * m
  
  if (d > 2) {
    lossD1 <- -colSums(kr(kr1, X, weight))
  } else if (d == 2) {
    lossD1 <- -colSums(kr1 * weight * X)
  }
  lossD1 <- matrix(lossD1, p, (d - 1))
  
  return(list(loss, lossD1))
}


setClassUnion("listOrMatrix", c("list", "matrix"))
setClass("LFIRlasso", representation(call = "call", data = "list", coefficients = "listOrMatrix", 
                                     logL = "numeric", BIC = "numeric", AIC = "numeric", Dof = "numeric", iter = "numeric", 
                                     maxlambda = "numeric", lambda = "numeric", distribution = "character"))
setClass("LFIRtune", representation(call = "call", select = "LFIRlasso", 
                                    path = "data.frame", select.list = "list"))

