#beta
LMNIRpen1 <- function(X, Y, base = NULL, penalty=NULL, trace = FALSE, n.factors=2, maxit = 1000) {
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
  out.list$n.factor <- n.factors
  out.list$iter=iter-1
  out.list$model.coefs$alpha <- new.model.coefs[1:p]
  out.list$model.coefs$phi <- matrix(new.model.coefs[(p+1):(p+k2*p)], k2, p)
  out.list$model.coefs$Beta <- matrix(new.model.coefs[-(1:((k2+1)*p))], k, p) 
  
  if(!is.null(penalty)){
    
    V.new <- cbind(V,f.hat)
    V.new <- t(t(V.new)/sqrt(colSums(V.new^2)))*sqrt(n)
    
    # sweep <- LFIRtune(formula=X~1+V.new, dist="MN", va.sigma=new.va.sigma, penalty=penalty, ngridpt=20, init=matrix(new.model.coefs, nrow = k+k2+1, ncol = p, byrow = T), penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    # sweep <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty=penalty, va.sigma=new.va.sigma, ngridpt=30, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    sweep1 <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty="lasso", va.sigma=0, ngridpt=30, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    sweep2 <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty="SCAD", va.sigma=0, ngridpt=30, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    sweep3 <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty="MCP", va.sigma=0, ngridpt=30, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    sweep4 <- LFIRtune(formula=X~1+V.new, dist="MN", k=k, penalty="group", va.sigma=0, ngridpt=60, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    # sweep <- LFIRtune(formula=X~1+V.new, dist="MN", penalty=penalty, va.sigma=new.va.sigma, lambdas=seq(from=300,to=350,by=1), penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    # sweep <- LFIRlasso(formula=X~1+V.new, dist="MN", k=k, penalty=penalty, lambda=0, va.sigma=new.va.sigma, penidx=c(FALSE,rep(TRUE,ncol(V.new))))
    # sweep4 <- MGLMtune(formula=X~1+V.new, dist="MN", penalty="sweep", ngridpt=30, penidx=c(FALSE,rep(TRUE,ncol(V.new))))

    new.model.coefs1 <- c(t(sweep1@select@coefficients))
    new.model.coefs2 <- c(t(sweep2@select@coefficients))
    new.model.coefs3 <- c(t(sweep3@select@coefficients))
    new.model.coefs4 <- c(t(sweep4@select@coefficients))
    # new.model.coefs <- c(t(sweep@coefficients))
    out.list$model.coefs1$alpha <- new.model.coefs1[1:p]
    out.list$model.coefs1$phi <- matrix(new.model.coefs1[(p+1):(p+k2*p)], k2, p ,byrow = T)
    out.list$model.coefs1$Beta <- matrix(new.model.coefs1[-(1:((k2+1)*p))], k, p,byrow = T)
    out.list$model.coefs2$alpha <- new.model.coefs2[1:p]
    out.list$model.coefs2$phi <- matrix(new.model.coefs2[(p+1):(p+k2*p)], k2, p,byrow = T)
    out.list$model.coefs2$Beta <- matrix(new.model.coefs2[-(1:((k2+1)*p))], k, p,byrow = T)
    out.list$model.coefs3$alpha <- new.model.coefs3[1:p]
    out.list$model.coefs3$phi <- matrix(new.model.coefs3[(p+1):(p+k2*p)], k2, p,byrow = T)
    out.list$model.coefs3$Beta <- matrix(new.model.coefs3[-(1:((k2+1)*p))], k, p,byrow = T)
    out.list$model.coefs4$alpha <- new.model.coefs4[1:p]
    out.list$model.coefs4$phi <- matrix(new.model.coefs4[(p+1):(p+k2*p)], k2, p,byrow = T)
    out.list$model.coefs4$Beta <- matrix(new.model.coefs4[-(1:((k2+1)*p))], k, p,byrow = T)
    out.list$lassosolution <- sweep1
    out.list$SCADsolution <- sweep2
    out.list$MCPsolution <- sweep3
    out.list$groupsolution <- sweep4
    # out.list$MGLMsolution <- sweep4
    
  }
  
  ## print the output
  
  out.list$va.coefs$Mu <- ff.hat
  # out.list$va.coefs$Sigma <- matrix(new.va.sigma^2, n, k)
  
  
  return(out.list)
}