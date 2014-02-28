##
## MCMC spatial kriging algorithm
##
## John Tipton
##
## Created 01.04.2014
## Last updated 01.22.2014
##
##
## Model: y_t = H_t %*% X %*% B_t + espilon_t
##
##        B_t ~ N(mu_B, Sigma_B) 
##
##        mu_B ~ N(mu_0, Sigma_0)  
##
## X_t = f(Y_tau) where the Y's are the simulated curves for years tau where the simulated data exits and f() is a PCA transformation
##
## y_t are the sampled observations for year t
##
## H_t is the selection matrix that ties observation locations for the sampled data to the simulated data
##
## mu_0 and Sigma_0 are hyperparameters
##

mcmc.1d <- function(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, s.star){
  
  ##
  ## Libraries and Subroutines
  ##
  
  ##
  ## Ideas to speed up code <- update cH.list after updating c and feed this into functions
  ## update Sig.mat <- tcH.list %*% C.star.inv %*% cH.list after updating c and C.star.inv and feed into functions
  ## store diagonal matrices as vectors and expand into matrices as needed
  ##
  
  ## changing C.star to C.star.inv and C.star.inv to C.star
  ##
  ## Predictive Process
  ##
  
  make.sum.sigma.beta <- function(s, beta, mu.beta, Lambda.inv){
   t(beta[, s] - mu.beta) %*% Lambda.inv %*% (beta[, s] - mu.beta)
  }
  
  ###############################
  ## This needs to be reviewed ##
  ###############################
  make.mh <- function(s, beta, Sigma.epsilon.tilde.plus.epsilon, Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.inv, C.star, C.star.inv, cH.list, tcH.list){
    ( - 1 / 2) * (determinant(C.star + cH.list[[s]] %*% Sigma.epsilon.tilde.plus.epsilon.inv[[s]] %*% tcH.list[[s]], logarithm = TRUE)$modulus[1] + determinant(C.star.inv, logarithm = TRUE)$modulus[1] + determinant(Sigma.epsilon.tilde.plus.epsilon[[s]], logarithm = TRUE)$modulus[1]) - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% (Sigma.inv[[s]]) %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
  }
  
  make.identity.list <- function(s, nt){
    if(length(nt) == 1){
      diag(nt)
    } else {
      diag(nt[s])
    }
  }
  
  
  make.Sigma.epsilon.tilde.plus.epsilon <- function(s, sigma.squared.epsilon, sigma.squared.eta, cH.list, tcH.list, C.star.inv, I.nt){
    (rep(sigma.squared.epsilon + sigma.squared.eta, nt[s]) - diag(tcH.list[[s]] %*% C.star.inv %*% cH.list[[s]])) * I.nt[[s]] 
  }
  
  ##
  ## Could speed this up if needed
  ##
  make.Sigma.epsilon.tilde.plus.epsilon.inv <- function(s, Sigma.epsilon.tilde.plus.epsilon){
    diag(1 / diag(Sigma.epsilon.tilde.plus.epsilon[[s]]))
  }
  
  make.c <- function(sigma.squared.eta, phi){
    sigma.squared.eta * exp( - D.0 / phi)
  }
  
  make.cH.list <- function(s, c, H.list){
    c[, H.list[[s]]]
  }
  
  make.tcH.list <- function(s, cH.list){
    t(cH.list[[s]])
  }
  
  make.C.star <- function(sigma.squared.eta, phi){
    sigma.squared.eta * exp( - D.star / phi)
  }
  
  make.C.star.inv <- function(C.star){
    solve(C.star)
  }
  
  make.Sigma <- function(s, Sigma.epsilon.tilde.plus.epsilon, cH.list, tcH.list, C.star.inv){
    tcH.list[[s]] %*% C.star.inv %*% cH.list[[s]] + Sigma.epsilon.tilde.plus.epsilon[[s]]
  }  
  
  make.Sigma.inv <- function(s, Sigma.epsilon.tilde.plus.epsilon.inv, C.star, cH.list, tcH.list){
    Sigma.epsilon.tilde.plus.epsilon.inv[[s]] - Sigma.epsilon.tilde.plus.epsilon.inv[[s]] %*% tcH.list[[s]] %*% solve(C.star + cH.list[[s]] %*% Sigma.epsilon.tilde.plus.epsilon.inv[[s]] %*% tcH.list[[s]]) %*% cH.list[[s]] %*% Sigma.epsilon.tilde.plus.epsilon.inv[[s]]
  }
  
  make.Sigma.beta.det <- function(s, sigma.squared.beta, tau){
    (1 / sigma.squared.beta)^tau
  }
  
  make.Sigma.full <- function(c, C.star.inv, sigma.squared.eta, sigma.squared.epsilon){
    Sig.mat <- t(c) %*% C.star.inv %*% c
    Sig.mat + (rep((sigma.squared.eta + sigma.squared.epsilon), ncells) - diag(Sig.mat)) * I.full
  }
  
  make.fort.batch <- function(s, beta, Sigma.full){
    Sigma.chol <- chol(Sigma.full)
    devs <- rnorm(ncells)
    if(dim(beta)[1] == 1){
      X * beta[s] + Sigma.chol %*% devs
    } else {
      X %*% beta[, s] + Sigma.chol %*% devs
    }
  }
  
  ##
  ## Initialize parameters
  ## 
  num.pca <- 3
  X.pca <- prcomp(X, center = FALSE, retx = TRUE)
  X <- X.pca$x[, 1:num.pca]  
  t <- length(Y.list)
  if(is.null(dim(X)) == TRUE){ncells <- length(X)} else {ncells <- dim(X)[1]}
  if(is.null(dim(X)) == TRUE){tau <- 1} else {tau <- dim(X)[2]}
  nt <- c() 
  for(s in 1:t){
    nt[s] <- length(Y.list[[s]])
  }
  nt.sum <- sum(nt)
  I.full <- diag(ncells)
  
  ## Distance Matrices for Predictive Process
  n.knots <- length(s.star) # predictive process
  D.star <- as.matrix(dist(s.star))
  D.0 <- matrix(nrow = n.knots, ncol = ncells)
  for(i in 1:n.knots){
    for(j in 1:ncells){
      D.0[i, j] <- sqrt((s.star[i] - locs[j])^2)
    }
  }
  
  ## Initialze process
  beta <- matrix(0, nrow = tau, ncol = t)
  
  ## Initialize parameter model
  lambda <- X.pca$sdev^2
  Lambda <- diag(lambda[1:num.pca])
  Lambda.determinant <- det(Lambda)
  Lambda.inv <- solve(Lambda)
  sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * Lambda
  Sigma.beta.inv <- 1 / sigma.squared.beta * Lambda
  
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  D <- as.matrix(dist(locs))
  
  vec.0 <- rep(0, n.knots)
  C.star <- make.C.star(sigma.squared.eta = sigma.squared.eta, phi = phi)
  eta.star <- rmvnorm(1, vec.0, C.star)
  C.star.inv <- make.C.star.inv(C.star)
  c <- make.c(sigma.squared.eta, phi)
  cH.list <- lapply(1:t, make.cH.list, c = c, H.list = H.list)
  tcH.list <- lapply(1:t, make.tcH.list, cH.list = cH.list)
  
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)
  Sigma.epsilon.tilde.plus.epsilon <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta, cH.list = cH.list, tcH.list = tcH.list, C.star.inv = C.star.inv, I.nt = I.nt)
  Sigma.epsilon.tilde.plus.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon)
  
  Sigma.0 <- sigma.squared.0 * diag(tau)
  Sigma.0.inv <- solve(Sigma.0)
  
  Sigma <- lapply(1:t, make.Sigma, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon, cH.list = cH.list, tcH.list = tcH.list, C.star.inv = C.star.inv) 
  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.inv, C.star = C.star, cH.list = cH.list, tcH.list = tcH.list)
  devs <- rnorm(tau)
  Sigma.chol <- chol(Sigma.0)
  mu.beta <- backsolve(Sigma.chol, backsolve(Sigma.0, mu.0, transpose = TRUE) + devs)
  
  n.burn <- floor(n.mcmc / 5) + 1
  fort.raster.batch <- matrix(0, ncells, t)   
  
  tHX.list <- vector('list', length = t) # speeds up computation by not calculating each MCMC step
  HX.list <- vector('list', length = t)
  for(s in 1:t){
    HX.list[[s]] <-   X[H.list[[s]], ]
    tHX.list[[s]] <- t(HX.list[[s]])
  }
  
  ##
  ## Initialize Storage
  ##
  
  beta.save <- array(dim = c(tau, t, n.mcmc))
  sigma.squared.beta.save <- vector(length = n.mcmc)
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  sigma.squared.eta.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc) 
  mu.beta.save <- matrix(NA, nrow = tau, ncol = n.mcmc)
  fort.raster <- matrix(0, nrow = ncells, ncol = t)
  var.save <- matrix(0, nrow = ncells, ncol = t)
  var.save.temp <- array(dim = c(100, ncells, t))
  phi.accept <- 0  
  eta.accept <- 0
  epsilon.accept <- 0
  
  ##
  ## Begin MCMC loop
  ##
  
  for(k in 1:n.mcmc){
    if(k %% 100 == 0) cat(" ", k) 
    
    ##
    ## Sample Beta
    ##
    
    for(s in 1:t){
      devs <- rnorm(tau)
      beta.A.chol <- chol(tHX.list[[s]] %*% Sigma.inv[[s]] %*% HX.list[[s]] + Sigma.beta.inv)
      beta.b <- tHX.list[[s]] %*% Sigma.inv[[s]] %*% Y.list[[s]] + Sigma.beta.inv %*% mu.beta
      beta[, s] <- backsolve(beta.A.chol, backsolve(beta.A.chol, beta.b, transpose = TRUE) + devs)
    }
    
    ##
    ## Sample mu.beta
    ##
    
    devs <- rnorm(tau)
    mu.beta.A.chol <- chol(t * Sigma.beta.inv + Sigma.0.inv)
    mu.beta.b <- apply(Sigma.beta.inv %*% beta, 1, sum) + Sigma.0.inv %*% mu.0
    mu.beta <- backsolve(mu.beta.A.chol, backsolve(mu.beta.A.chol, mu.beta.b, transpose = TRUE) + devs)
    
    ##
    ## Sample sigma.squared.beta
    ##
    
    sigma.squared.beta <- 1 / rgamma(1, alpha.beta + t * tau / 2, beta.beta + 1 / 2 * sum(sapply(1:t, make.sum.sigma.beta, beta = beta, mu.beta = mu.beta, Lambda.inv = Lambda.inv)))
    Sigma.beta <- sigma.squared.beta * Lambda
    Sigma.beta.inv <- 1 / sigma.squared.beta * Lambda.inv
    
    ##
    ## Sample sigma.squared.eta
    ##
    
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      c.star <- make.c(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      cH.list.star <- lapply(1:t, make.cH.list, c = c.star, H.list = H.list)
      tcH.list.star <- lapply(1:t, make.tcH.list, cH.list.star)
      C.star.star <- make.C.star(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      C.star.star.inv <- make.C.star.inv(C.star = C.star.star)
      Sigma.epsilon.tilde.plus.epsilon.star <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta.star, cH.list = cH.list.star, tcH.list = tcH.list.star, C.star.inv = C.star.star.inv, I.nt = I.nt)
      Sigma.epsilon.tilde.plus.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, cH.list = cH.list.star, tcH.list = tcH.list.star, C.star.inv = C.star.star.inv)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, C.star = C.star.star, cH.list = cH.list.star, tcH.list = tcH.list.star)
      mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, Sigma.inv = Sigma.star.inv, C.star = C.star.star, C.star.inv = C.star.star.inv, cH.list = cH.list.star, tcH.list = tcH.list.star)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
      mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, cH.list = cH.list, tcH.list = tcH.list)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
      mh.eta <- exp(mh.eta.1 - mh.eta.2)
      
      if(mh.eta > runif(1)){
        sigma.squared.eta <- sigma.squared.eta.star
        c <- c.star
        cH.list <- cH.list.star
        tcH.list <- tcH.list.star
        C.star <- C.star.star
        C.star.inv <- C.star.star.inv
        Sigma.epsilon.tilde.plus.epsilon <- Sigma.epsilon.tilde.plus.epsilon.star
        Sigma.epsilon.tilde.plus.epsilon.inv <- Sigma.epsilon.tilde.plus.epsilon.star.inv
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        eta.accept <- eta.accept + 1 / n.mcmc
      }
      rm(c.star)
      rm(cH.list.star)
      rm(tcH.list.star)
      rm(C.star.star)
      rm(C.star.star.inv)
      rm(Sigma.epsilon.tilde.plus.epsilon.star)
      rm(Sigma.epsilon.tilde.plus.epsilon.star.inv)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(sigma.squared.eta.star)
    
    ##
    ## Sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon.star <- rnorm(1, sigma.squared.epsilon, sigma.squared.epsilon.tune)
    if(sigma.squared.epsilon.star > 0){
      Sigma.epsilon.tilde.plus.epsilon.star <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon, sigma.squared.epsilon = sigma.squared.epsilon.star, sigma.squared.eta = sigma.squared.eta, cH.list = cH.list, tcH.list = tcH.list, C.star.inv = C.star.inv, I.nt = I.nt)
      Sigma.epsilon.tilde.plus.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, cH.list = cH.list, tcH.list = tcH.list, C.star.inv = C.star.inv)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, C.star = C.star, cH.list = cH.list, tcH.list = tcH.list)
      mh.epsilon.1 <-  sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, Sigma.inv = Sigma.star.inv, C.star = C.star, C.star.inv = C.star.inv, cH.list = cH.list, tcH.list = tcH.list)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, cH.list = cH.list, tcH.list = tcH.list)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)
      
      if(mh.epsilon > runif(1)){
        sigma.squared.epsilon <- sigma.squared.epsilon.star
        Sigma.epsilon.tilde.plus.epsilon <- Sigma.epsilon.tilde.plus.epsilon.star
        Sigma.epsilon.tilde.plus.epsilon.inv <- Sigma.epsilon.tilde.plus.epsilon.star.inv
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        epsilon.accept <- epsilon.accept + 1 / n.mcmc
      }
      rm(Sigma.epsilon.tilde.plus.epsilon.star)
      rm(Sigma.epsilon.tilde.plus.epsilon.star.inv)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(sigma.squared.epsilon.star)
    
    ##
    ## Sample phi
    ##
    
    phi.star <- rnorm(1, phi, phi.tune)
    if(phi.star > 0){
      c.star <- make.c(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
      cH.list.star <- lapply(1:t, make.cH.list, c = c.star, H.list = H.list)
      tcH.list.star <- lapply(1:t, make.tcH.list, cH.list = cH.list.star)
      C.star.star <- make.C.star(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
      C.star.star.inv <- make.C.star.inv(C.star = C.star.star)
      Sigma.epsilon.tilde.plus.epsilon.star <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta, cH.list = cH.list.star, tcH.list = tcH.list.star, C.star.inv = C.star.star.inv, I.nt = I.nt)
      Sigma.epsilon.tilde.plus.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, cH.list = cH.list.star, tcH.list = tcH.list.star, C.star.inv = C.star.star.inv)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, C.star = C.star.star, cH.list = cH.list.star, tcH.list = tcH.list.star)        
      mh.phi.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon.star, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.star.inv, Sigma.inv = Sigma.star.inv, C.star = C.star.star, C.star.inv = C.star.star.inv, cH.list = cH.list.star, tcH.list = tcH.list.star)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
      mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.tilde.plus.epsilon = Sigma.epsilon.tilde.plus.epsilon, Sigma.epsilon.tilde.plus.epsilon.inv = Sigma.epsilon.tilde.plus.epsilon.inv, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, cH.list = cH.list, tcH.list = tcH.list)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
      
      if(mh.phi > runif(1)){
        phi <- phi.star
        c <- c.star
        cH.list <- cH.list.star
        tcH.list <- tcH.list.star
        C.star <- C.star.star
        C.star.inv <- C.star.star.inv
        Sigma.epsilon.tilde.plus.epsilon <- Sigma.epsilon.tilde.plus.epsilon.star
        Sigma.epsilon.tilde.plus.epsilon.inv <- Sigma.epsilon.tilde.plus.epsilon.star.inv
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        phi.accept <- phi.accept + 1 / n.mcmc
      }
      rm(c.star)
      rm(cH.list.star)
      rm(tcH.list.star)
      rm(C.star.star)
      rm(C.star.star.inv)
      rm(Sigma.epsilon.tilde.plus.epsilon.star)
      rm(Sigma.epsilon.tilde.plus.epsilon.star.inv)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(phi.star)
    
    ##
    ## Simulate random field
    ##
    
    if(k > n.burn){
      if(k %% 10 == 0){
        Sigma.full <- make.Sigma.full(c, C.star.inv, sigma.squared.eta, sigma.squared.epsilon)
        fort.raster <- fort.raster + 10 / (n.mcmc - n.burn) * sapply(1:t, make.fort.batch, beta = beta, Sigma.full = Sigma.full)
        if(k %% 1000 == 0){
          var.save.temp[100, , ] <- fort.raster
        } else {
          var.save.temp[(k %% 1000) / 10, , ] <- fort.raster
        }
      }
      if(k %% 1000 == 0){
        for(s in 1:t){
          var.save[, s] <- var.save[, s] + apply(var.save.temp[, , s], 2, var) / ((n.mcmc - n.burn) / 1000)
          var.save.temp <- array(0, dim = c(100, ncells, t))
        }
      }
    }
    
    ##
    ## Save variables
    ## 
    
    beta.save[, , k] <- beta
    sigma.squared.beta.save[k] <- sigma.squared.beta
    sigma.squared.eta.save[k] <- sigma.squared.eta
    sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
    mu.beta.save[, k] <- mu.beta  
    phi.save[k] <- phi	
  }
  
  ##
  ## Write output
  ##
  
  list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, sigma.squared.eta.save = sigma.squared.eta.save, mu.beta.save = mu.beta.save, n.mcmc = n.mcmc, fort.raster = fort.raster, phi.accept = phi.accept, eta.accept = eta.accept, epsilon.accept = epsilon.accept, phi.save = phi.save, var.save = var.save)#, MSPE.save = MSPE.save)
}
