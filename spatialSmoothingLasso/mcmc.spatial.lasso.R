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

## params = c(mu.0, Sigma.0, alpha.beta, beta.beta, alpha.eta, beta.eta, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, n.mcmc = 5000)
  
  ##
  ## Libraries and Subroutines
  ##

  #source('dinvgamma.R')
	
  make.sum.sigma.beta <- function(s, beta, mu.beta){
    t(beta[, s] - mu.beta) %*% (beta[, s] - mu.beta)
  }
    
  make.D.list <- function(s, H.list, locs){
    as.matrix(dist(locs[H.list[[s]]]))
  }
  
  make.R.list <- function(s, sigma.squared.eta, phi, D.list){
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
  
  #make.Sigma.epsilon.inv <- function(s, sigma.squared.epsilon, I.nt){
  #  1 / sigma.squared.epsilon * I.nt[[s]]
  #}
    
  make.Sigma <- function(s, R.list, Sigma.epsilon){
    R.list[[s]] + Sigma.epsilon[[s]]
  }
  
  make.Sigma.inv <- function(s, Sigma){
    solve(Sigma[[s]])
  }
  
	make.mh <- function(s, beta, Sigma, Sigma.inv){
	  ( - 1 / 2) * determinant(Sigma[[s]], logarithm = TRUE)$modulus[1] - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% Sigma.inv[[s]] %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
	}

  make.fort.batch <- function(s, beta, H.list, Y.list, Sigma.full, ncells){
    temp <- vector(length = ncells)
    temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + Sigma.full[ - H.list[[s]], H.list[[s]]] %*% solve(Sigma.full[H.list[[s]], H.list[[s]]]) %*% (Y.list[[s]] - X[H.list[[s]], ] %*% beta[, s])
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
  
  ## Initialize parameter model
  sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  I.beta <- diag(tau)
  Sigma.beta <- sigma.squared.beta * I.beta
  Sigma.beta.inv <- solve(Sigma.beta)
  
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  D <- as.matrix(dist(locs))
  D.list <- lapply(1:t, make.D.list, H.list = H.list, locs = locs)
  R.list <- lapply(1:t, make.R.list, sigma.squared.eta = sigma.squared.eta, phi = phi, D.list = D.list)
    
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)
  I.full <- diag(ncells)
  Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)
  #Sigma.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)

  Sigma.0 <- sigma.squared.0 * I.beta
  Sigma.0.inv <- solve(Sigma.0)

  Sigma <- lapply(1:t, make.Sigma, R.list = R.list, Sigma.epsilon = Sigma.epsilon)
  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma)
  
  devs <- rnorm(t)
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
  	
    sigma.squared.beta <- 1 / rgamma(1, alpha.beta + nt.sum / 2, beta.beta + 1 / 2 * sum(sapply(1:t, make.sum.sigma.beta, beta = beta, mu.beta = mu.beta)))
  	Sigma.beta <- sigma.squared.beta * I.beta
  	Sigma.beta.inv <- solve(Sigma.beta)
  	
  	##
    ## Sample sigma.squared.eta
    ##
  	
  	sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      R.list.star <- lapply(1:t, make.R.list, sigma.squared.eta = sigma.squared.eta.star, phi = phi, D.list = D.list)
      Sigma.star <- lapply(1:t, make.Sigma, R.list = R.list.star, Sigma.epsilon = Sigma.epsilon)
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
      #Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, sigma.squared.epsilon = sigma.squared.epsilon, I.nt = I.nt)
      Sigma.star <- lapply(1:t, make.Sigma, R.list = R.list, Sigma.epsilon = Sigma.epsilon.star)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.epsilon.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)
  	  	
  	  if(mh.epsilon > runif(1)){
        sigma.squared.epsilon <- sigma.squared.epsilon.star
  	    Sigma.epsilon <- Sigma.epsilon.star
  	  	#Sigma.epsilon.inv <- Sigma.epsilon.star.inv
  	  	Sigma <- Sigma.star
  	    Sigma.inv <- Sigma.star.inv
        epsilon.accept <- epsilon.accept + 1 / n.mcmc
  	  }
  	  rm(Sigma.epsilon.star)
  	  #rm(Sigma.epsilon.star.inv)
    	rm(Sigma.star)
    	rm(Sigma.star.inv)
  	}
    rm(sigma.squared.epsilon.star)
  	
  	##
  	## Sample phi
  	##
  	
    phi.star <- rnorm(1, phi, phi.tune)
  	if(phi.star > 0){
      R.list.star <- lapply(1:t, make.R.list, sigma.squared.eta = sigma.squared.eta, phi = phi.star, D.list = D.list)	  
  	  Sigma.star <- lapply(1:t, make.Sigma, R.list = R.list.star, Sigma.epsilon = Sigma.epsilon)
  	  Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma = Sigma.star)
      mh.phi.1 <-	sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma.star, Sigma.inv = Sigma.star.inv)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
  		mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma = Sigma, Sigma.inv = Sigma.inv)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
  		mh.phi <- exp(mh.phi.1 - mh.phi.2)
  		  		
  	if(mh.phi > runif(1)){
      phi <- phi.star
      R.list <- R.list.star
      Sigma <- Sigma.star
      Sigma.inv <- Sigma.star.inv
      phi.accept <- phi.accept + 1 / n.mcmc
  	}
    rm(R.list.star)
    rm(Sigma.star)
    rm(Sigma.star.inv)
  }
  rm(phi.star) 
  	
  ##
  ## Simulate random field
  ##
    
  if(k > n.burn){
    if(k %% 100 == 0){
      for(s in 1:t){
        var.save[, s] <- var.save[, s] + apply(var.save.temp[, , s], 2, var) / ((n.mcmc - n.burn) / 100)
        fort.raster.batch <- matrix(0, ncells, t)
      }
      var.save.temp <- array(0, dim = c(100, ncells, t))
    }
#     Sigma.full <- (sigma.squared.eta * exp( - D / phi))# + sigma.squared.epsilon * diag(ncells))
    Sigma.full <- (sigma.squared.eta * exp( - D / phi)) + sigma.squared.epsilon * I.full
    fort.raster <- fort.raster + 1 / (n.mcmc - n.burn) * sapply(1:t, make.fort.batch, beta = beta, H.list = H.list, Y.list = Y.list, Sigma.full = Sigma.full, ncells = ncells)
    #chol.Sigma <- chol(Sigma.full)
    #fort.raster <- fort.raster + (1 / (n.mcmc - n.burn)) * sapply(1:t, make.fort.batch, beta = beta, chol.Sigma = chol.Sigma, ncells = ncells)
    #rmvnorm(t(X) %*% beta, Sigma.full)
    var.save.temp[k %% 100 + 1, , ] <- fort.raster
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
