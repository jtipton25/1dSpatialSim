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

mcmc.1d <- function(Y.list, H.list, X, locs, n.mcmc, mu.0, Sigma.0, alpha.epsilon, beta.epsilon, alpha.beta, beta.beta, alpha.phi, beta.phi, mu.beta, sigma.squared.eta.tune, sigma.squared.epsilon.tune, phi.tune){
  
  
  ##
  ## Libraries and Subroutines
  ##
  
  #source('dinvgamma.R')
  
  make.sum.sigma.beta <- function(s, beta, mu.beta){
    t(beta[, s] - mu.beta) %*% Lambda.inv %*% (beta[, s] - mu.beta)
  }
  
  make.D.list <- function(s, H.list, locs){
    as.matrix(dist(locs[H.list[[s]]]))
  }
  
  make.Sigma.eta <- function(sigma.squared.eta, phi){
    #   make.R.list <- function(s, sigma.squared.eta, phi, D.list, LtL.list){
    #     sigma.squared.eta * LtL.list[[s]] %*% exp( - D.list[[s]] / phi) %*% LtL.list[[s]]
    #     sigma.squared.eta * tL.list[[s]] %*% exp( - D.list[[s]] / phi) %*% L.list[[s]]
    sigma.squared.eta * t(L) %*% exp( - D / phi) %*% L
  }
  
  #   make.L.list <- function(s, I.nt, HX.list, tHX.list){
  #     e <- eigen(I.nt[[s]] - HX.list[[s]] %*% 
  #     solve(tHX.list[[s]] %*% HX.list[[s]])
  #      %*% tHX.list[[s]])
  #     idx <- round(e$values, 4) == 1
  #     return(e$vectors[, idx])
  #   }
  #   
  #   make.tL.list <- function(s, L.list){
  #     return(t(L.list[[s]]))
  #   }
  # 
  #    make.LtL.list <- function(s, L.list, tL.list){
  #      return(L.list[[s]] %*% tL.list[[s]])
  #    }
  #    
  #    make.tLL.list <- function(s, L.list, tL.list){
  #      return(tL.list[[s]] %*% L.list[[s]])
  #    }
  #   
  make.identity.list <- function(s, nt){
    if(length(nt) == 1){
      diag(nt)
    } else {
      diag(nt[s])
    }
  }
  
  make.Sigma.epsilon <- function(s, sigma.squared.epsilon){
    sigma.squared.epsilon * I.nt[[s]]
  }
  
  make.Sigma.epsilon.inv <- function(s, sigma.squared.epsilon){
    1 / sigma.squared.epsilon * I.nt[[s]]
  }
  
  make.Sigma <- function(Sigma.eta, sigma.squared.epsilon){
    L %*% Sigma.eta %*% tL + sigma.squared.epsilon * I.full
  }
  
  make.Sigma.inv <- function(Sigma){
    solve(Sigma)
  }
  
  make.eta <- function(s, Sigma.eta){
    #    rmvnorm(1, sigma = Sigma.eta[[s]])
    chol(Sigma.eta) %*% rnorm(n.L)
  }
  
  make.mh <- function(s, beta, eta, sigma.squared.epsilon){
    # 	  ( - 1 / 2) * determinant(Sigma.epsilon[[s]], logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s] - L.list[[s]] %*% eta[[s]]) %*% Sigma.epsilon.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s] - L.list[[s]] %*% eta[[s]])
    ( - 1 / 2) * sigma.squared.epsilon^(- nt[s] / 2) * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s] - (L %*% eta[[s]])[H.list[[s]]]) %*%  (Y.list[[s]] - HX.list[[s]] %*% beta[, s] - (L %*% eta[[s]])[H.list[[s]]]) / sigma.squared.epsilon
  }
  
  #   make.fort.batch <- function(s, beta, Sigma.full, Sigma.inv){
  #     temp <- vector(length = ncells)
  #     temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + Sigma.full[ - H.list[[s]], H.list[[s]]] %*% Sigma.inv[H.list[[s]], H.list[[s]]] %*% (Y.list[[s]] - X[H.list[[s]], ] %*% beta[, s])
  #     temp[H.list[[s]]] <- Y.list[[s]]
  #     return(temp)
  #   }
  
  #   make.fort.batch <- function(s, beta, Sigma.full, Sigma.inv){
  #     temp <- vector(length = ncells)
  #     temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + Sigma.full[ - H.list[[s]], H.list[[s]]] %*% Sigma.inv[H.list[[s]], H.list[[s]]] %*% (Y.list[[s]] - X[H.list[[s]], ] %*% beta[, s])
  #     temp[H.list[[s]]] <- Y.list[[s]]
  #     return(temp)
  #   }
  make.fort.batch <- function(s, beta, eta, sigma.squared.epsilon){
    devs <- rnorm(ncells - nt[s], mean = 0, sd = sqrt(sigma.squared.epsilon))
    temp <- vector(length = ncells)
    temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + (L %*% eta[[s]])[ - H.list[[s]]] + devs
    temp[H.list[[s]]] <- Y.list[[s]]
    return(temp)
  }
  
  ##
  ## Initialize parameters
  ## 
  
  t <- length(Y.list)
  if(is.null(dim(X)) == TRUE){ncells <- length(X)} else {ncells <- dim(X)[1]}
  if(is.null(dim(X)) == TRUE){tau <- 1} else {tau <- dim(X)[2]}
  e <- eigen(diag(ncells) - X %*% solve(t(X) %*% X) %*% t(X))
  idx <- round(e$values, 4) == 1
  L <- e$vectors[, idx]
  tL <- t(L)
  tLL <- t(L) %*% L
  n.L <- dim(L)[2]
  nt <- c() 
  for(s in 1:t){
    nt[s] <- length(Y.list[[s]])
  }
  nt.sum <- sum(nt)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)
  I.tau <- diag(tau)
  I.full <- diag(ncells)
  tHX.list <- vector('list', length = t) # speeds up computation by not calculating each MCMC step
  HX.list <- vector('list', length = t)
  for(s in 1:t){
    HX.list[[s]] <-   X[H.list[[s]], ]
    tHX.list[[s]] <- t(HX.list[[s]])
  }
  #   L.list <- lapply(1:t, make.L.list, I.nt = I.nt, HX.list = HX.list, tHX.list = tHX.list)
  #   tL.list <- lapply(1:t, make.tL.list, L.list = L.list)
  #   LtL.list <- lapply(1:t, make.LtL.list, L.list = L.list, tL.list = tL.list)
  #   tLL.list <- lapply(1:t, make.tLL.list, L.list = L.list, tL.list = tL.list)  
  #   ## Initialze process
  beta <- matrix(0, nrow = tau, ncol = t)
  
  ## Initialize parameter model
  sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * I.tau
  Sigma.beta.inv <- 1 / sigma.squared.beta * I.tau
  
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  D <- as.matrix(dist(locs))
  D.list <- lapply(1:t, make.D.list, H.list = H.list, locs = locs)
  #   Sigma.eta <- lapply(1:t, make.Sigma.eta, sigma.squared.eta = sigma.squared.eta, phi = phi)
  Sigma.eta <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta, phi = phi)
  #   n.L.list <- vector(length = t)
  #   for(s in 1:t){
  #     n.L.list[s] <- dim(Sigma.eta[[s]])[1]  
  #   }
  
  eta <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta)
  eta.star <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta)
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon)
  Sigma.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon)
  
  #  Sigma.0 <- sigma.squared.0 * diag(tau)
  Sigma.0.inv <- solve(Sigma.0)
  
  #  Sigma <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta, Sigma.epsilon = Sigma.epsilon)
  #  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma)
  
  ##
  ## fix this
  ##
  
  devs <- rnorm(tau)
  Sigma.chol <- chol(Sigma.0)
  mu.beta <- backsolve(Sigma.chol, backsolve(Sigma.chol, mu.0, transpose = TRUE) + devs)
  
  n.burn <- floor(n.mcmc / 5) + 1
  fort.raster.batch <- matrix(0, ncells, t)   
  
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
      beta.A.chol <- chol(tHX.list[[s]] %*% Sigma.epsilon.inv[[s]] %*% HX.list[[s]] + Sigma.beta.inv)
      beta.b <- tHX.list[[s]] %*% Sigma.epsilon.inv[[s]] %*% (Y.list[[s]] - (L %*% eta[[s]])[H.list[[s]]]) + Sigma.beta.inv %*% mu.beta
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
    
    sigma.squared.beta <- 1 / rgamma(1, alpha.beta + t * tau / 2, beta.beta + 1 / 2 * sum(sapply(1:t, make.sum.sigma.beta, beta = beta, mu.beta = mu.beta)))
    Sigma.beta <- sigma.squared.beta * Lambda
    Sigma.beta.inv <- 1 / sigma.squared.beta * Lambda.inv
    
    ##
    ## Sample eta
    ##
    
    A.chol <- chol(tLL / sigma.squared.epsilon + Sigma.eta)
    for(s in 1:t){
      devs <- rnorm(n.L)
      b <- 1 / sigma.squared.epsilon * tL[, H.list[[s]]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
      eta[[s]] <- backsolve(A.chol, backsolve(A.chol, b, transpose = TRUE) + devs)
    }
    
    
    ##
    ## Sample sigma.squared.eta
    ##
    
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      #       Sigma.eta.star <- lapply(1:t, make.Sigma.eta, sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      Sigma.eta.star <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      A.chol <- chol(tLL / sigma.squared.epsilon + Sigma.eta.star)
      for(s in 1:t){
        devs <- rnorm(n.L)
        b <- 1 / sigma.squared.epsilon * tL[, H.list[[s]]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
        eta.star[[s]] <- backsolve(A.chol, backsolve(A.chol, b, transpose = TRUE) + devs)
      }
      #       Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta.star, Sigma.epsilon = Sigma.epsilon)
      #       Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star) 
      mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta.star, sigma.squared.epsilon = sigma.squared.epsilon)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
      mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, sigma.squared.epsilon = sigma.squared.epsilon)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
      mh <- exp(mh.eta.1 - mh.eta.2)
      
      if(mh > runif(1)){
        sigma.squared.eta <- sigma.squared.eta.star
        Sigma.eta <- Sigma.eta.star
        eta <- eta.star
        #         Sigma <- Sigma.star
        #         Sigma.inv <- Sigma.star.inv
        eta.accept <- eta.accept + 1 / n.mcmc
      }
      #       rm(Sigma.star)
      #       rm(Sigma.star.inv)
    }
    rm(sigma.squared.eta.star)
    
    ##
    ## Sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon.star <- rnorm(1, sigma.squared.epsilon, sigma.squared.epsilon.tune)
    if(sigma.squared.epsilon.star > 0){
      #       Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon.star)
      #       Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon.star)
      #       Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta, Sigma.epsilon = Sigma.epsilon.star)
      #       Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.epsilon.1 <-	sum(sapply(1:t, make.mh, beta = beta, eta = eta, sigma.squared.epsilon = sigma.squared.epsilon.star)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, sigma.squared.epsilon = sigma.squared.epsilon.star)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)
      
      if(mh.epsilon > runif(1)){
        sigma.squared.epsilon <- sigma.squared.epsilon.star
        #   	    Sigma.epsilon <- Sigma.epsilon.star
        #Sigma.epsilon.inv <- Sigma.epsilon.star.inv
        #   	  	Sigma <- Sigma.star
        #   	    Sigma.inv <- Sigma.star.inv
        epsilon.accept <- epsilon.accept + 1 / n.mcmc
      }
      #   	  rm(Sigma.epsilon.star)
      #rm(Sigma.epsilon.star.inv)
      #     	rm(Sigma.star)
      #     	rm(Sigma.star.inv)
    }
    rm(sigma.squared.epsilon.star)
    
    ##
    ## Sample phi
    ##
    
    phi.star <- rnorm(1, phi, phi.tune)
    if(phi.star > 0){
      Sigma.eta.star <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
      A.chol <- chol(1 / sigma.squared.epsilon * tLL + Sigma.eta.star)
      for(s in 1:t){
        devs <- rnorm(n.L)
        b <- 1 / sigma.squared.epsilon * tL[, H.list[[s]]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
        eta.star[[s]] <- backsolve(A.chol, backsolve(A.chol, b, transpose = TRUE) + devs)
      }
      #       Sigma.star <- lapply(1:t, make.Sigma, Sigma.eta = Sigma.eta.star, Sigma.epsilon = Sigma.epsilon)
      #   	  Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.phi.1 <-	sum(sapply(1:t, make.mh, beta = beta, eta = eta.star, sigma.squared.epsilon = sigma.squared.epsilon)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
      mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, sigma.squared.epsilon = sigma.squared.epsilon)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
      
      if(mh.phi > runif(1)){
        phi <- phi.star
        Sigma.eta <- Sigma.eta.star
        eta <- eta.star
        #       Sigma <- Sigma.star
        #       Sigma.inv <- Sigma.star.inv
        phi.accept <- phi.accept + 1 / n.mcmc
      }
      rm(Sigma.eta.star)
      #     rm(Sigma.star)
      #     rm(Sigma.star.inv)
    }
    rm(phi.star) 
    
    ##
    ## Simulate random field
    ##
    
    ##
    ## Update this section of code
    ##
    
    if(k > n.burn){
      if(k %% 10 == 0){
        ##
        ## check this
        ##
            Sigma.full <- (sigma.squared.eta * exp( - D / phi)) + sigma.squared.epsilon * I.full
        #     Sigma <- make.Sigma(Sigma.eta = Sigma.eta, sigma.squared.epsilon = sigma.squared.epsilon)
        #     Sigma.inv <- make.Sigma.inv(Sigma = Sigma)
        #     fort.raster.tmp <- sapply(1:t, make.fort.batch, beta = beta, Sigma.full = Sigma.full, Sigma.inv = Sigma.inv)
        fort.raster.tmp <- sapply(1:t, make.fort.batch, beta = beta, eta = eta, sigma.squared.epsilon = sigma.squared.epsilon)
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
