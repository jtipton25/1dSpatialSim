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
    t(beta[, s] - mu.beta) %*% (beta[, s] - mu.beta)
  }
  
#   make.D.list <- function(s, H.list, locs){
#     as.matrix(dist(locs[H.list[[s]]]))
#   }
  
#   make.Sigma.eta <- function(sigma.squared.eta, phi){
# #     sigma.squared.eta * tL %*% exp( - D / phi) %*% L
#     sigma.squared.eta * tL %*% exp( - D / phi) %*% L
#   }
#   
#   make.eta <- function(s, Sigma.eta){
# #    rmvnorm(1, sigma = Sigma.eta[[s]])
#     chol(Sigma.eta) %*% rnorm(n.L)
#   }

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

  make.Sigma.fixed <- function(phi.vec){
    n.phi <- length(phi.vec)
    tmp <- vector('list', length = n.phi)
    for(i in 1:n.phi){
      tmp[[i]] <- LtL %*% exp( - D / phi) %*% LtL
    }
  return(tmp)
  }
  
  make.Sigma.full <- function(sigma.squared.eta, sigma.squared.epsilon, phi){
#     sigma.squared.eta * LtL %*% exp( - D / phi) %*% LtL + sigma.squared.epsilon * I.full
    idx.sig <- which(phi == phi.vec)
    sigma.squared.eta * Sigma.fixed[[idx.sig]] + sigma.squared.epsilon * I.full
  }

  make.Sigma <- function(s, Sigma.full){
    Sigma.full[H.list[[s]], H.list[[s]]]
#     Sigma.eta + sigma.squared.epsilon * I.full
  }
  
  make.Sigma.inv <- function(s, Sigma){
    solve(Sigma[[s]])
  }
  
#   make.eta <- function(s, Sigma.eta){
#     chol(Sigma.eta) %*% rnorm(n.L)
#   }
  
#   make.mh <- function(s, beta, eta, Sigma.epsilon, Sigma.epsilon.inv){
  make.mh <- function(s, beta, Sigma, Sigma.inv){
#     	  ( - 1 / 2) * determinant(Sigma, logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s] - L[H.list[[s]], ] %*% eta) %*% Sigma.inv[H.list[[s]], H.list[[s]]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
#   ( - 1 / 2) * determinant(Sigma.epsilon[[s]], logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s] - (L %*% eta[[s]])[H.list[[s]]]) %*% Sigma.epsilon.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s] - (L %*% eta[[s]])[H.list[[s]]])
    ( - 1 / 2) * determinant(Sigma[[s]], logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% Sigma.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
  }
  
#     make.fort.batch <- function(s, beta, sigma.squared.epsilon){
#       X %*% beta[, s] + L %*% eta[[s]] + rnorm(ncells, 0, sqrt(sigma.squared.epsilon))
#     }

  make.fort.batch <- function(s, beta, Sigma.full, Sigma.inv){
    tmp <- vector(length = ncells)
    tmp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + Sigma.full[ - H.list[[s]], H.list[[s]]] %*% Sigma.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
    tmp[H.list[[s]]] <- Y.list[[s]]
    return(tmp)
  }

  ##
  ## Initialize parameters
  ## 
  phi.vec <- (1:50) / 50
  t <- length(Y.list)
  if(is.null(dim(X)) == TRUE){ncells <- length(X)} else {ncells <- dim(X)[1]}
  if(is.null(dim(X)) == TRUE){tau <- 1} else {tau <- dim(X)[2]}
  e <- eigen(diag(ncells) - X %*% solve(t(X) %*% X) %*% t(X))
  idx <- round(e$values, 4) == 1
  L <- e$vectors[, idx]
  tL <- t(L)
#   tLL <- tL %*% L
  LtL <- L %*% tL
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

  beta <- matrix(0, nrow = tau, ncol = t)
  
  ## Initialize parameter model
  sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * I.tau
  Sigma.beta.inv <- 1 / sigma.squared.beta * I.tau
  
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
#   phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  phi <- sample(phi.vec, 1)
  D <- as.matrix(dist(locs))
#   D.list <- lapply(1:t, make.D.list, H.list = H.list, locs = locs)
#   Sigma.eta <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta, phi = phi)
  
#   eta <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta)
#   eta.star <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta)
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  Sigma.fixed <- make.Sigma.fixed(phi.vec)
  Sigma.full <- make.Sigma.full(sigma.squared.eta = sigma.squared.eta, sigma.squared.epsilon = sigma.squared.epsilon, phi = phi)
  Sigma <- lapply(1:t, make.Sigma, Sigma.full = Sigma.full)
  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma)
#   Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon)
#   Sigma.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon)
  
  #  Sigma.0 <- sigma.squared.0 * diag(tau)
  Sigma.0.inv <- solve(Sigma.0)
  
#   Sigma <- make.Sigma(Sigma.eta = Sigma.eta, sigma.squared.epsilon = sigma.squared.epsilon)
#   Sigma.inv <- make.Sigma.inv(Sigma = Sigma)
  
  ##
  ## fix this
  ##
  
  devs <- rnorm(tau)
  Sigma.chol <- chol(Sigma.0)
  mu.beta <- backsolve(Sigma.chol, backsolve(Sigma.chol, mu.0, transpose = TRUE) + devs)
  
  n.burn <- floor(n.mcmc / 5) + 1
 
  
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
  fort.raster.batch <- matrix(0, ncells, t)  
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
    
    sigma.squared.beta <- 1 / rgamma(1, alpha.beta + t * tau / 2, beta.beta + 1 / 2 * sum(sapply(1:t, make.sum.sigma.beta, beta = beta, mu.beta = mu.beta)))
    Sigma.beta <- sigma.squared.beta * I.tau
    Sigma.beta.inv <- 1 / sigma.squared.beta * I.tau
    
    ##
    ## Sample sigma.squared.eta
    ##
    
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      Sigma.full.star <- make.Sigma.full(sigma.squared.eta = sigma.squared.eta.star, sigma.squared.epsilon = sigma.squared.epsilon, phi = phi)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.full = Sigma.full.star)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
#       Sigma.eta.star <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
#       eta.star <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta.star)
      #Sigma.star <- make.Sigma(Sigma.eta = Sigma.eta.star, sigma.squared.epsilon = sigma.squared.epsilon)
      #Sigma.star.inv <- make.Sigma.inv(Sigma = Sigma.star) 
        mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta)
#       mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta.star, Sigma.epsilon = Sigma.epsilon, Sigma.epsilon.inv = Sigma.epsilon.inv)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
#       mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
        mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta)
#       mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, Sigma.epsilon = Sigma.epsilon, Sigma.epsilon.inv = Sigma.epsilon.inv)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
#       mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
      mh <- exp(mh.eta.1 - mh.eta.2)
      
      if(mh > runif(1)){
        sigma.squared.eta <- sigma.squared.eta.star
#         Sigma.eta <- Sigma.eta.star
        Sigma.full <- Sigma.full.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        eta.accept <- eta.accept + 1 / n.mcmc
      }
#       rm(Sigma.eta.star)
      rm(Sigma.full.star)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(sigma.squared.eta.star)
    
    ##
    ## Sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon.star <- rnorm(1, sigma.squared.epsilon, sigma.squared.epsilon.tune)
    if(sigma.squared.epsilon.star > 0){
      Sigma.full.star <- make.Sigma.full(sigma.squared.eta = sigma.squared.eta, sigma.squared.epsilon = sigma.squared.epsilon.star, phi = phi)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.full = Sigma.full.star)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
#       Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon.star)
#       Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon.star)
#       Sigma.star <- make.Sigma(Sigma.eta = Sigma.eta, sigma.squared.epsilon = sigma.squared.epsilon)
      mh.epsilon.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon)
#       mh.epsilon.1 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, Sigma.epsilon = Sigma.epsilon.star, Sigma.epsilon.inv = Sigma.epsilon.star.inv)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
#       Sigma.star.inv <- make.Sigma.inv(Sigma = Sigma.star)
#       mh.epsilon.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon)
# mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, Sigma.epsilon = Sigma.epsilon, Sigma.epsilon.inv = Sigma.epsilon.inv)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
#       mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)
      
      if(mh.epsilon > runif(1)){
        sigma.squared.epsilon <- sigma.squared.epsilon.star
#         Sigma.epsilon <- Sigma.epsilon.star
#         Sigma.epsilon.inv <- Sigma.epsilon.star.inv
        Sigma.full <- Sigma.full.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        epsilon.accept <- epsilon.accept + 1 / n.mcmc
      }
      rm(Sigma.full.star)
      rm(Sigma.star)
      rm(Sigma.star.inv)
    }
    rm(sigma.squared.epsilon.star)
    
    ##
    ## Sample phi
    ##
    
#     phi.star <- rnorm(1, phi, phi.tune)
    phi.idx <- which(phi == phi.vec)
    proposal <- sample( ( - phi.tune:phi.tune)[ - (phi.tune + 1)], 1)

    if(phi.idx + proposal <= length(phi.vec) && phi.idx + proposal > 0){
      phi.star <- phi.vec[phi.idx + proposal]
      Sigma.full.star <- make.Sigma.full(sigma.squared.eta = sigma.squared.eta, sigma.squared.epsilon = sigma.squared.epsilon, phi = phi.star)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.full = Sigma.full.star)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.phi.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) #+ dinvgamma(phi.star, alpha.phi, beta.phi)
      mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) #+ dinvgamma(phi, alpha.phi, beta.phi)
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
      
      if(mh.phi > runif(1)){
        phi <- phi.star
        Sigma.full <- Sigma.full.star
        Sigma <- Sigma.star
        Sigma.inv <- Sigma.star.inv
        phi.accept <- phi.accept + 1 / n.mcmc
      }
      #       rm(Sigma.eta.star)
      rm(Sigma.full.star)
      rm(Sigma.star)
      rm(Sigma.star.inv)
      rm(phi.star) 
    } 

#     if(phi.star > 0){
#       Sigma.full.star <- make.Sigma.full(sigma.squared.eta = sigma.squared.eta, sigma.squared.epsilon = sigma.squared.epsilon, phi = phi.star)
#       Sigma.star <- lapply(1:t, make.Sigma, Sigma.full = Sigma.full.star)
#       Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
# #       Sigma.eta.star <- make.Sigma.eta(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
# #       eta.star <- lapply(1:t, make.eta, Sigma.eta = Sigma.eta.star)
# # #       Sigma.star <- make.Sigma(Sigma.eta = Sigma.eta.star, sigma.squared.epsilon = sigma.squared.epsilon)
# #       Sigma.star.inv <- make.Sigma.inv(Sigma = Sigma.star)
#       mh.phi.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(phi.star, alpha.phi, beta.phi)
# #       mh.phi.1 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta.star, Sigma.epsilon = Sigma.epsilon, Sigma.epsilon.inv = Sigma.epsilon.inv)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
# #       mh.phi.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
#       mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(phi, alpha.phi, beta.phi)
# #       mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, eta = eta, Sigma.epsilon = Sigma.epsilon, Sigma.epsilon.inv = Sigma.epsilon.inv)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
# #       mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
#       mh.phi <- exp(mh.phi.1 - mh.phi.2)
#       
#       if(mh.phi > runif(1)){
#         phi <- phi.star
# #         Sigma.eta <- Sigma.eta.star
#         Sigma.full <- Sigma.full.star
#         Sigma <- Sigma.star
#         Sigma.inv <- Sigma.star.inv
#         phi.accept <- phi.accept + 1 / n.mcmc
#       }
# #       rm(Sigma.eta.star)
#       rm(Sigma.full.star)
#       rm(Sigma.star)
#       rm(Sigma.star.inv)
#     }
#     rm(phi.star) 
#     
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
#         Sigma.full <- LtL %*% (sigma.squared.eta * exp( - D / phi)) %*% LtL + sigma.squared.epsilon * I.full
#         Sigma.full.inv <- solve(Sigma.full)
        fort.raster.tmp <- sapply(1:t, make.fort.batch, beta = beta, Sigma.full = Sigma.full, Sigma.inv = Sigma.inv)
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
