##
## MCMC principal component spatial kriging algorithm
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

  make.sum.sigma.beta <- function(s, beta, mu.beta){
    t(beta[, s] - mu.beta) %*% Lambda.inv %*% (beta[, s] - mu.beta)
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
  
  make.Sigma.epsilon.inv <- function(s, sigma.squared.epsilon, I.nt){
   1 / sigma.squared.epsilon * I.nt[[s]]
  }
    
  make.sum.sigma.epsilon <- function(s, beta){
    t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
  }

  ##
  ## Initialize parameters
  ## 

#   num.pca <- 3
  X.pca <- prcomp(X, center = TRUE, retx = TRUE)
  X <- X.pca$x
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
  
  ## Initialize parameter model
  lambda <- X.pca$sdev^2
  Lambda <- diag(lambda)
  Lambda.determinant <- det(Lambda)
  Lambda.inv <- solve(Lambda)
  U <- X.pca$rotation
  H.hat <- U %*% (Lambda)^(1 / 2)
  sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * Lambda
  Sigma.beta.inv <- 1 / sigma.squared.beta * Lambda.inv
  
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)
  I.full <- diag(ncells)
  Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)
  Sigma.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)

  
  Sigma.0.inv <- solve(Sigma.0)

  mu.beta <- backsolve(chol(Sigma.0), backsolve(chol(Sigma.0), mu.0, transpose = TRUE) + rnorm(tau))

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
  mu.beta.save <- matrix(NA, nrow = tau, ncol = n.mcmc)
  fort.raster <- matrix(0, nrow = ncells, ncol = t)
  var.save <- matrix(0, nrow = ncells, ncol = t)
  var.save.temp <- array(dim = c(100, ncells, t))
  
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
      beta.b <- tHX.list[[s]] %*% Sigma.epsilon.inv[[s]] %*% Y.list[[s]] + Sigma.beta.inv %*% mu.beta
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
    ## Sample sigma.squared.epsilon
    ##
  	
    sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon + nt.sum, beta.epsilon + 1 / 2 * sum(sapply(1:t, make.sum.sigma.epsilon, beta = beta)))

    ##
    ## Simulate random field
    ##
    
    if(k > n.burn){
      if(k %% 10 == 0){
        fort.raster.tmp <- X %*% beta
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
  sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
  mu.beta.save[, k] <- mu.beta  
  }
  
##
## Write output
##
  
list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, 
     sigma.squared.epsilon.save = sigma.squared.epsilon.save, mu.beta.save = mu.beta.save, 
     fort.raster = fort.raster, var.save = var.save)#, MSPE.save = MSPE.save)
}
