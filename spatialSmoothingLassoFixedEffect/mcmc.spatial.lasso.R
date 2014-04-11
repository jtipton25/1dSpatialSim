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

mcmc.1d <- function(Y.list, H.list, X, locs, n.mcmc, mu.beta, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, alpha.lambda, beta.lambda, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune, gamma.squared.tune){
    
  ##
  ## Libraries and Subroutines
  ##
  
  make.sum.beta <- function(s, beta){
    t(beta[, s] - mu.beta) %*% (beta[, s] - mu.beta)
  }
   
  make.D.list <- function(s, H.list, locs){
    as.matrix(dist(locs[H.list[[s]]]))
  }
  
  make.Sigma.eta <- function(s, sigma.squared.eta, phi, D.list){
    sigma.squared.eta * exp( - D.list[[s]] / phi)
  }
  
  make.identity.list <- function(s, nt){
    if(length(nt) == 1){
      diag(nt)
    } else {
      diag(nt[s])
    }
  }
  
  make.Sigma.epsilon <- function(s, sigma.squared.epsilon, I.nt){
    sigma.squared.epsilon * I.nt[[s]]
  }
  
  make.Sigma <- function(s, Sigma.eta, Sigma.epsilon){
    Sigma.eta[[s]] + Sigma.epsilon[[s]]
  }
  
  make.Sigma.inv <- function(s, Sigma){
    solve(Sigma[[s]])
  }
  
  make.mh <- function(s, beta, Sigma, Sigma.inv){
    ( - 1 / 2) * determinant(Sigma[[s]], logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% Sigma.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
  }
  
  make.mh.epsilon <- function(s, beta, Sigma, Sigma.inv, Sigma.beta, Sigma.beta.inv){
    ( - 1 / 2) * determinant(Sigma[[s]],  logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% Sigma.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s]) + ( - 1 / 2) * determinant(Sigma.beta,  logarithm = TRUE)$modulus[1] - 1 / 2 * beta[, s] %*% Sigma.beta.inv %*% beta[, s]
  }
  
#   make.mh.gamma <- function(s, beta, sigma.squared.epsilon, gamma.squared, lambda.squared){
#     ( - T / 2) * log(sigma.squared.epsilon * gamma.squared[s]) - sum(sapply(1:t, make.sum.beta, beta = beta)) / (2 * sigma.squared.epsilon + gamma.squared[s]) - lambda.squared * gamma.squared[s] / 2
#   }

  make.mh.gamma <- function(s, beta, gamma.squared, Sigma.beta.inv){
    ( - 1 / 2) * log(prod(gamma.squared)) + ( - 1 / 2) * beta[, s] %*% Sigma.beta.inv %*% beta[, s]
  }
  
  make.fort.batch <- function(s, beta, c.Y, Sigma){
    temp <- vector(length = ncells)
    temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + c.Y[ - H.list[[s]], H.list[[s]]] %*% solve(Sigma[[s]]) %*% (Y.list[[s]] - X[H.list[[s]], ] %*% beta[, s])
    temp[H.list[[s]]] <- Y.list[[s]]
    return(temp)
  }
  
  ##
  ## Initialize parameters
  ## 
  
  t <- length(Y.list)
  if(is.null(dim(X)) == TRUE){ncells <- length(X)} else {ncells <- dim(X)[1]}
  if(is.null(dim(X)) == TRUE){tau <- 1} else {tau <- dim(X)[2]}
  nt <- c() 
  for(s in 1:t){
    nt[s] <- length(Y.list[[s]])
  }
  nt.sum <- sum(nt)
  
  ## Initialze process
  beta <- matrix(0, nrow = tau, ncol = t)
  
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  D <- as.matrix(dist(locs))
  D.list <- lapply(1:t, make.D.list, H.list = H.list, locs = locs)
  Sigma.eta <- lapply(1:t, make.Sigma.eta, sigma.squared.eta = sigma.squared.eta, phi = phi, D.list = D.list)
  
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)
  I.full <- diag(ncells)
  Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)

  Sigma <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta, Sigma.epsilon = Sigma.epsilon)
  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma)

  lambda.squared <- 2 * rgamma(1, alpha.lambda, beta.lambda)
  gamma.squared <- rexp(tau, lambda.squared / 2)
  D.gamma <- diag(gamma.squared)
  D.gamma.inv <- diag(1 / gamma.squared)
  Sigma.beta <- sigma.squared.epsilon * D.gamma
  Sigma.beta.inv <- 1 / sigma.squared.epsilon * D.gamma.inv

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
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  sigma.squared.eta.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc) 
  gamma.squared.save <- matrix(NA, nrow = tau, ncol = n.mcmc)
  lambda.squared.save <- vector(length = n.mcmc)
  fort.raster <- matrix(0, nrow = ncells, ncol = t)
  var.save <- matrix(0, nrow = ncells, ncol = t)
  var.save.temp <- array(dim = c(100, ncells, t))
  phi.accept <- 0  
  eta.accept <- 0
  epsilon.accept <- 0
  gamma.accept <- 0
  
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
    ## Sample sigma.squared.eta
    ##
    
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      Sigma.eta.star <- lapply(1:t, make.Sigma.eta, sigma.squared.eta = sigma.squared.eta.star, phi = phi, D.list = D.list)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta.star, Sigma.epsilon = Sigma.epsilon)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
      mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
      mh <- exp(mh.eta.1 - mh.eta.2)
      
      if(mh > runif(1)){
        sigma.squared.eta <- sigma.squared.eta.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        eta.accept <- eta.accept + 1 / n.mcmc
      }
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(sigma.squared.eta.star)
    
    ##
    ## Sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon.star <- rnorm(1, sigma.squared.epsilon, sigma.squared.epsilon.tune)
    if(sigma.squared.epsilon.star > 0){
      Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon.star, I.nt = I.nt)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta, Sigma.epsilon = Sigma.epsilon.star)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      Sigma.beta.star <- sigma.squared.epsilon.star * D.gamma
      Sigma.beta.star.inv <- 1 / sigma.squared.epsilon.star * D.gamma.inv
      mh.epsilon.1 <-	sum(sapply(1:t, make.mh.epsilon, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv, Sigma.beta  = Sigma.beta.star, Sigma.beta.inv = Sigma.beta.star.inv)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh.epsilon, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv, Sigma.beta  = Sigma.beta, Sigma.beta.inv = Sigma.beta.inv)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)
      
      if(mh.epsilon > runif(1)){
        sigma.squared.epsilon <- sigma.squared.epsilon.star
        Sigma.epsilon <- Sigma.epsilon.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        Sigma.beta <- Sigma.beta.star
        Sigma.beta.inv <- Sigma.beta.inv
        epsilon.accept <- epsilon.accept + 1 / n.mcmc
      }
      rm(Sigma.epsilon.star)
      rm(Sigma.star)
      rm(Sigma.star.inv)
      rm(Sigma.beta.star)
      rm(Sigma.beta.star.inv)
    }
    rm(sigma.squared.epsilon.star)
    
    ##
    ## Sample phi
    ##
    
    phi.star <- rnorm(1, phi, phi.tune)
    if(phi.star > 0){
      Sigma.eta.star <- lapply(1:t, make.Sigma.eta, sigma.squared.eta = sigma.squared.eta, phi = phi.star, D.list = D.list)	  
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta.star, Sigma.epsilon = Sigma.epsilon)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.phi.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
      mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
      
      if(mh.phi > runif(1)){
        phi <- phi.star
        Sigma.eta <- Sigma.eta.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        phi.accept <- phi.accept + 1 / n.mcmc
      }
      rm(Sigma.eta.star)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(phi.star) 
    
    ##
    ## Sample gamma.squared
    ##

##
## double check this
##
    gamma.squared.star <- rnorm(tau, gamma.squared, rep(sigma.squared.gamma.tune, tau))
    if(min(gamma.squared.star) > 0){
      D.gamma.star <- diag(gamma.squared.star)
      D.gamma.star.inv <- diag(1 / gamma.squared.star) 
      Sigma.beta.star <- sigma.squared.epsilon * D.gamma.star
      Sigma.beta.star.inv <- 1 / sigma.squared.epsilon * D.gamma.star.inv
      mh.gamma.1 <- sum(sapply(1:t, make.mh.gamma, beta = beta, gamma.squared = gamma.squared.star, Sigma.beta.inv = Sigma.beta.star.inv)) - sum((lambda.squared / 2) * gamma.squared.star)
      mh.gamma.2 <- sum(sapply(1:t, make.mh.gamma, beta = beta, gamma.squared = gamma.squared, Sigma.beta.inv = Sigma.beta.inv)) - sum((lambda.squared / 2) * gamma.squared)
      mh.gamma <- exp(mh.gamma.1 - mh.gamma.2)
    
      if(mh.gamma > runif(1)){
        gamma.squared <- gamma.squared.star
        D.gamma <- D.gamma.star
        D.gamma.inv <- D.gamma.star.inv
        Sigma.beta <- Sigma.beta.star
        Sigma.beta.inv <- Sigma.beta.star.inv
        gamma.accept <- gamma.accept + 1 / n.mcmc
      }
    }
    rm(gamma.squared.star)

    ##
    ## sample lambda.squared
    ##
  
    lambda.squared <- rgamma(1, alpha.lambda + tau, beta.lambda / 2 + sum(gamma.squared) / 2)
    
    ##
    ## Simulate random field
    ##
    
    if(k > n.burn){
      if(k %% 10 == 0){
        c.Y <- sigma.squared.eta * exp( - D / phi)
        fort.raster.tmp <- sapply(1:t, make.fort.batch, beta = beta, c.Y = c.Y, Sigma = Sigma)
        fort.raster <- fort.raster + 10 / (n.mcmc - n.burn) * fort.raster.tmp
        if(k %% 1000 == 0){
          var.save.temp[100, , ] <- fort.raster.tmp
        } else {
          var.save.temp[(k %% 1000) / 10, , ] <- fort.raster.tmp
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
#     sigma.squared.beta.save[k] <- sigma.squared.beta
    sigma.squared.eta.save[k] <- sigma.squared.eta
    sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
    lambda.squared.save[k] <- lambda.squared
    gamma.squared.save[, k] <- gamma.squared
    phi.save[k] <- phi	
  }
  
  ##
  ## Write output
  ##
  
  list(beta.save = beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, sigma.squared.eta.save = sigma.squared.eta.save, lambda.squared.save = lambda.squared.save, n.mcmc = n.mcmc, fort.raster = fort.raster, phi.accept = phi.accept, eta.accept = eta.accept, epsilon.accept = epsilon.accept, phi.save = phi.save, gamma.accept = gamma.accept, gamma.squared.save = gamma.squared.save, lambda.squared.save = lambda.squared.save, var.save = var.save)#, MSPE.save = MSPE.save)
}
