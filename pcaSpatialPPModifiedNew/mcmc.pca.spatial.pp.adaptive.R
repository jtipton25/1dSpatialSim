##
## MCMC principal component spatial predictive process algorithm
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

  make.mh <- function(s, beta, Sigma.epsilon.inv, Sigma.epsilon, Sigma.inv, C.star, C.star.inv, c){
    cH.list <- c[, H.list[[s]]]
    ( - 1 / 2) * (determinant(C.star + cH.list %*% Sigma.epsilon.inv[[s]] %*% t(cH.list), logarithm = TRUE)$modulus[1] + determinant(C.star.inv, logarithm = TRUE)$modulus[1] + determinant(Sigma.epsilon[[s]])$modulus[1]) - 1 / 2 * t(Y.list[[s]] - HX.list[[s]] %*% beta[, s]) %*% (Sigma.inv[[s]]) %*% (Y.list[[s]] - HX.list[[s]] %*% beta[, s])
  }
    
    make.identity.list <- function(s, nt){
      if(length(nt) == 1){
        diag(nt)
      } else {
        diag(nt[s])
      }
    }
    
   make.sum.sigma.beta <- function(s, beta, mu.beta, Lambda.inv){
     t(beta[, s] - mu.beta) %*% Lambda.inv %*% (beta[, s] - mu.beta)
   }
  
  make.Sigma.epsilon <- function(s, sigma.squared.epsilon, sigma.squared.eta, c, C.star.inv, I.nt){
    cH.list <- c[, H.list[[s]]]
    (sigma.squared.epsilon + sigma.squared.eta - diag((t(cH.list) %*% C.star.inv %*% cH.list))) * I.nt[[s]]
  }
    
  make.Sigma.epsilon.inv <- function(s, Sigma.epsilon){
    solve(Sigma.epsilon[[s]])
  }
  
  make.c <- function(sigma.squared.eta, phi){
    sigma.squared.eta * exp( - D.0 / phi)
  }
  
  make.C.star <- function(sigma.squared.eta, phi){
    sigma.squared.eta * exp( - D.star / phi)
  }
  
  make.C.star.inv <- function(C.star){
    solve(C.star)
  }
  
  make.Sigma <- function(s, Sigma.epsilon, C.star.inv, c, H.list){
    cH.list <- c[, H.list[[s]]]
    (t(cH.list) %*% C.star.inv %*% cH.list) + Sigma.epsilon[[s]]
  }
  
  make.Sigma.inv <- function(s, Sigma.epsilon.inv, C.star, c, H.list){
    cH.list <- c[, H.list[[s]]]
    Sigma.epsilon.inv[[s]] - Sigma.epsilon.inv[[s]] %*% t(cH.list) %*% solve(C.star + cH.list %*% Sigma.epsilon.inv[[s]] %*% t(cH.list)) %*% cH.list %*% Sigma.epsilon.inv[[s]]
  }
  
  make.Sigma.beta.det <- function(s, sigma.squared.beta, tau){
    (1 / sigma.squared.beta)^tau * Lambda.determinant
  }
  
#   make.fort.batch <- function(s, beta, c, C.star, C.star.inv, sigma.squared.epsilon){
#     w.star <- rmvnorm(1, vec.0, C.star) # w.star should be c.star using the notation in the writeup but c.star is taken by the mh step
#     w.tilde <- t(c) %*% C.star.inv %*% t(w.star) #likewise, w.tilde is c.tilde in the writeup
#   	if(dim(beta)[1] == 1){
#   		X * beta[s] + w.tilde + rnorm(ncells, 0, sigma.squared.epsilon)
#   	} else {
#   		X %*% beta[, s] + w.tilde + rnorm(ncells, 0, sigma.squared.epsilon)
#   	}
#   }

  make.fort.batch <- function(s, beta, Sigma.inv, Sigma.full){
    temp <- vector(length = ncells)
    temp[ - H.list[[s]]] <- X[ - H.list[[s]], ] %*% beta[, s] + Sigma.full[ - H.list[[s]], H.list[[s]]] %*% Sigma.inv[[s]] %*% (Y.list[[s]] - X[H.list[[s]], ] %*% beta[, s])
    temp[H.list[[s]]] <- Y.list[[s]]  
    return(temp)
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

  sigma.squared.beta <- 0.00998
#   sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * Lambda
  Sigma.beta.inv <- 1 / sigma.squared.beta * Lambda.inv
  
  sigma.squared.eta <- 3.679
  phi <- 1.041
#   sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
#   phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  D <- as.matrix(dist(locs))
    
  sigma.squared.epsilon <- 0.0569
#   sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  I.nt <- lapply(1:t, make.identity.list, nt = nt)

  Sigma.0.inv <- solve(Sigma.0)

  vec.0 <- rep(0, n.knots)
  C.star <- make.C.star(sigma.squared.eta = sigma.squared.eta, phi = phi)
  w.star <- rmvnorm(1, vec.0, C.star)
  C.star.inv <- make.C.star.inv(C.star)
  c <- make.c(sigma.squared.eta, phi)
  Sigma.epsilon <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta, c = c, C.star.inv = C.star.inv, I.nt = I.nt)
  Sigma.epsilon.inv <- lapply(1:t, make.Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon)

  Sigma <- lapply(1:t, make.Sigma, Sigma.epsilon = Sigma.epsilon, C.star.inv = C.star.inv, c = c, H.list = H.list) 
  Sigma.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.inv = Sigma.epsilon.inv, C.star = C.star, c = c, H.list = H.list)
  
#   devs <- rnorm(tau)
#   Sigma.chol <- chol(Sigma.0)
  mu.beta <- c(-0.159, 0.138, 0.00180)
#   mu.beta <- backsolve(Sigma.chol, backsolve(Sigma.0, mu.0, transpose = TRUE) + devs)

  n.burn <- floor(n.mcmc / 5) + 1
  fort.raster.batch <- matrix(0, ncells, t)   

  tHX.list <- vector('list', length = t) # speeds up computation by not calculating each MCMC step
  HX.list <- vector('list', length = t)
  for(s in 1:t){
    if(length(H.list[[s]]) == 1){
      HX.list[[s]] <- t(X[H.list[[s]], ])
      tHX.list[[s]] <- t(HX.list[[s]])
    } else {
      HX.list[[s]] <- X[H.list[[s]], ]
      tHX.list[[s]] <- t(HX.list[[s]])
    }
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
  eta.accept <- 0
  eta.accept.tmp <- 0
  epsilon.accept <- 0
  epsilon.accept.tmp <- 0
  phi.accept <- 0  
  phi.accept.tmp <- 0  

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
      beta.A.chol <-  chol(tHX.list[[s]] %*% Sigma.inv[[s]] %*% HX.list[[s]] + Sigma.beta.inv)
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

    if(k %% 100 == 0){
      if(eta.accept.tmp < 0.20){
        sigma.squared.eta.tune <- sigma.squared.eta.tune * 0.75
      } else if(eta.accept.tmp > 0.60){
        sigma.squared.eta.tune <- sigma.squared.eta.tune * 1.25
      }
      eta.accept.tmp <- 0
  	}
    
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.squared.eta.tune)
    if(sigma.squared.eta.star > 0){
      c.star <- make.c(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      C.star.star <- make.C.star(sigma.squared.eta = sigma.squared.eta.star, phi = phi)
      C.star.star.inv <- make.C.star.inv(C.star = C.star.star)
      Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta.star, c = c.star, C.star.inv = C.star.star.inv, I.nt = I.nt)
      Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon = Sigma.epsilon.star, C.star.inv = C.star.star.inv, c = c.star, H.list = H.list)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.inv = Sigma.epsilon.star.inv, C.star = C.star.star, c = c.star, H.list = H.list)
      mh.eta.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.star.inv, Sigma.epsilon = Sigma.epsilon.star, Sigma.inv = Sigma.star.inv, C.star = C.star.star, C.star.inv = C.star.star.inv, c = c.star)) + dinvgamma(sigma.squared.eta.star, alpha.eta, beta.eta, log = TRUE)
      mh.eta.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, c = c)) + dinvgamma(sigma.squared.eta, alpha.eta, beta.eta, log = TRUE)
      mh.eta <- exp(mh.eta.1 - mh.eta.2)
    	  	
       if(mh.eta > runif(1)){
         sigma.squared.eta <- sigma.squared.eta.star
         c <- c.star
         C.star <- C.star.star
         C.star.inv <- C.star.star.inv
         Sigma.epsilon <- Sigma.epsilon.star
         Sigma.epsilon.inv <- Sigma.epsilon.star.inv
         Sigma <- Sigma.star
         Sigma.inv <- Sigma.star.inv
         eta.accept <- eta.accept + 1 / n.mcmc
         eta.accept.tmp <- eta.accept.tmp + 1 / 100
        }
      }
    
    ##
    ## Sample sigma.squared.epsilon
    ##
  	
    if(k %% 100 == 0){
      if(epsilon.accept.tmp < 0.20){
        sigma.squared.epsilon.tune <- sigma.squared.epsilon.tune * 0.75
      } else if(eta.accept.tmp > 0.60){
        sigma.squared.epsilon.tune <- sigma.squared.epsilon.tune * 1.25
      }
      epsilon.accept.tmp <- 0
    }
    
    sigma.squared.epsilon.star <- rnorm(1, sigma.squared.epsilon, sigma.squared.epsilon.tune)
    if(sigma.squared.epsilon.star > 0){
      Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon.star, sigma.squared.eta = sigma.squared.eta, c = c, C.star.inv = C.star.inv, I.nt = I.nt)
      Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon)
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon = Sigma.epsilon.star, C.star.inv = C.star.inv, c = c, H.list = H.list)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv, Sigma.epsilon.inv = Sigma.epsilon.star.inv, C.star = C.star, c = c, H.list = H.list)
      mh.epsilon.1 <-  sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.star.inv, Sigma.epsilon = Sigma.epsilon.star, Sigma.inv = Sigma.star.inv, C.star = C.star, C.star.inv = C.star.inv, c = c)) + dinvgamma(sigma.squared.epsilon.star, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, c = c)) + dinvgamma(sigma.squared.epsilon, alpha.epsilon, beta.epsilon, log = TRUE)
      mh.epsilon <- exp(mh.epsilon.1 - mh.epsilon.2)

      if(mh.epsilon > runif(1)){
       	sigma.squared.epsilon <- sigma.squared.epsilon.star
       	Sigma.epsilon <- Sigma.epsilon.star
     	  Sigma.epsilon.inv <- Sigma.epsilon.star.inv
     	  Sigma <- Sigma.star
       	Sigma.inv <- Sigma.star.inv
       	epsilon.accept <- epsilon.accept + 1 / n.mcmc
       	epsilon.accept.tmp <- epsilon.accept.tmp + 1 / 100
      }
    }
    
    ##
    ## Sample phi
    ##
  	
    if(k %% 100 == 0){
      if(phi.accept.tmp < 0.20){
        phi.tune <- phi.tune * 0.75
      } else if(phi.accept.tmp > 0.60){
        phi.tune <- phi.tune * 1.25
      }
      phi.accept.tmp <- 0
    }
        
    phi.star <- rnorm(1, phi, phi.tune)
    if(phi.star > 0){
      c.star <- make.c(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
      C.star.star <- make.C.star(sigma.squared.eta = sigma.squared.eta, phi = phi.star)
      C.star.star.inv <- make.C.star.inv(C.star = C.star.star)
      Sigma.epsilon.star <- lapply(1:t, make.Sigma.epsilon, sigma.squared.epsilon = sigma.squared.epsilon, sigma.squared.eta = sigma.squared.eta, c = c.star, C.star.inv = C.star.star.inv, I.nt = I.nt)
      Sigma.epsilon.star.inv <- lapply(1:t, make.Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon)
      
      Sigma.star <- lapply(1:t, make.Sigma, Sigma.epsilon = Sigma.epsilon.star, C.star.inv = C.star.star.inv, c = c.star, H.list = H.list)
      Sigma.star.inv <- lapply(1:t, make.Sigma.inv,  Sigma.epsilon.inv = Sigma.epsilon.star.inv, C.star = C.star.star, c = c.star, H.list = H.list)
      mh.phi.1 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon, Sigma.inv = Sigma.star.inv, C.star = C.star.star, C.star.inv = C.star.star.inv, c = c.star)) + dinvgamma(phi.star, alpha.phi, beta.phi, log = TRUE)
      mh.phi.2 <- sum(sapply(1:t, make.mh, beta = beta, Sigma.epsilon.inv = Sigma.epsilon.inv, Sigma.epsilon = Sigma.epsilon, Sigma.inv = Sigma.inv, C.star = C.star, C.star.inv = C.star.inv, c = c)) + dinvgamma(phi, alpha.phi, beta.phi, log = TRUE)
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
     	  
      if(mh.phi > runif(1)){
        phi <- phi.star
     	  c <- c.star
     	  C.star <- C.star.star
     	  C.star.inv <- C.star.star.inv
     	  Sigma <- Sigma.star
       	Sigma.inv <- Sigma.star.inv
      	phi.accept <- phi.accept + 1 / n.mcmc
        phi.accept.tmp <- phi.accept.tmp + 1 / 100
      }
    }
    
    ##
    ## Simulate random field
    ##
    
    if(k > n.burn){
      if(k %% 10 == 0){
        Sigma.full <- sigma.squared.eta * exp( - D / phi)
        fort.raster.tmp <- sapply(1:t, make.fort.batch, beta = beta, Sigma.inv = Sigma.inv, Sigma.full = Sigma.full)
        fort.raster <- fort.raster + 1 / ((n.mcmc - n.burn) / 10) * fort.raster.tmp
        
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
## MSPE for cross-validation
##

#      temp <- c()
#      for(j in 1:t){
#        if(is.null(H.list[[j]]) == FALSE){
#          temp[j] <- (fort.raster[,j][H.list[[j]]] - y.val[,j][loc.id.val[[j]]])^2 
#        }
#      }
#      MSPE.save <- MSPE.save + mean(na.omit(temp)) / (n.mcmc - n.burn)
    
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
  
list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, sigma.squared.eta.save = sigma.squared.eta.save, mu.beta.save = mu.beta.save, n.mcmc = n.mcmc, fort.raster = fort.raster, phi.accept = phi.accept, eta.accept = eta.accept, epsilon.accept = epsilon.accept, phi.save = phi.save, var.save = var.save, sigma.squared.eta.tune = sigma.squared.eta.tune, sigma.squared.epsilon.tune = sigma.squared.epsilon.tune, phi.tune = phi.tune)#, MSPE.save = MSPE.save)
}

