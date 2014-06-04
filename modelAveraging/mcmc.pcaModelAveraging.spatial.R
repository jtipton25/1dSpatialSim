##
## MCMC sampler for orthogonal data augmentation
##

mcmc.pcaMA <- function(Y.list, X.o, H.list, params, tune, epsilon = 0.001){ #Y.new, X.new, for log scoring rule
  
  ##
  ## functions and subroutines
  ##
  
  make.mh <- function(i, sigma.squared, Sigma.t, Sigma.t.inv, gamma, beta.tilde.gamma){
    if(sum(gamma[[i]]) == 0){
      - n.o[i] / 2 * sigma.squared - 1 / 2 * determinant(Sigma.t[[i]], logarithm = TRUE)[1]$mod[1] - 1 / (2 * sigma.squared) * t(Y.list[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]]) 
    } else {
      - n.o[i] / 2 * sigma.squared - 1 / 2 * determinant(Sigma.t[[i]], logarithm = TRUE)[1]$mod[1] - 1 / (2 * sigma.squared) * t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]])  
    }
  }
  
  make.gamma.mh <- function(i, gamma, beta.hat, Sigma.full.inv, Y.c, tX.Sigma.full.inv.X, tX.Sigma.full.inv){
    sum(gamma[[i]] * log(pi.prior) + (1 - gamma[[i]]) * log(1 - pi.prior)) - 1 / (2 * sigma.squared) * (t(gamma[[i]] * beta.hat[[i]]) %*% tX.Sigma.full.inv.X %*% (gamma[[i]] * beta.hat[[i]]) - 2 * t(gamma[[i]] * beta.hat[[i]]) %*% tX.Sigma.full.inv %*% Y.c[[i]] + t(gamma[[i]] * beta.hat[[i]]) %*% Lambda %*% (gamma[[i]] * beta.hat[[i]]))
  }
  
  ##
  ## initialize fixed values
  ##
  
  n.mcmc <- params$n.mcmc
  alpha <- params$alpha
  pi.prior <- params$pi.prior
  lambda <- params$lambda
  alpha.eta <- params$alpha.eta
  beta.eta <- params$beta.eta
  phi.lower <- params$phi.lower
  phi.upper <- params$phi.upper
  D <- params$D
  
#   sigma.tune <- tune$sigma.tune
  phi.tune <- tune$phi.tune
  sigma.eta.tune <- tune$sigma.eta.tune
  gama.tune <- tune$gamma.tune
  
  t <- length(Y.list)
  X.pca <- prcomp(X.o)
  X <- X.pca$x
  tX <- t(X)
  delta <- X.pca$sdev^2
  p <- dim(X)[2]
  m <- dim(X)[1]
  n.o <- vector(length = t)
  n.u <- vector(length = t)
  for(i in 1:t){
    n.o[i] <- length(Y.list[[i]])
    n.u[i] <- m - n.o[i]
  }
  I.full <- diag(m)
  I.o <- vector('list', length = t)
  I.u <- vector('list', length = t)
  for(i in 1:t){
    I.o[[i]] <- diag(n.o[i])
    I.u[[i]] <- diag(n.u[i])
  }
  ## initialize random values
  
  ##
  ## choose a better starting value once the code is up and running
  ##
  
  sigma.squared <- 1
  
  ##
  Psi <- vector('list', length = t)
  gamma <- vector('list', length = t)
  for(i in 1:t){
    gamma[[i]] <- rbinom(p, 1, pi.prior)
  }
  gamma.star <- gamma
  Lambda <- diag(lambda)
  Lambda.gamma <- vector('list', length = t)
  for(i in 1:t){
    if(sum(gamma[[i]]) == 0){
      Lambda.gamma[[i]] <- 0
    } else {
  Lambda.gamma[[i]] <- diag(lambda[gamma[[i]] == 1])  
    }
  }
    
  ## 
  H.u.list <- vector('list', length = t)  
  for(i in 1:t){
    H.u.list[[i]] <- (1:m)[ - H.list[[i]]]
  }

  HX.o.list <- vector('list', length = t)
  tHX.o.list <- vector('list', length = t)
  HX.u.list <- vector('list', length = t)
  tHX.u.list <- vector('list', length = t)
  for(i in 1:t){
    HX.o.list[[i]] <- X[H.list[[i]], ]
    tHX.o.list[[i]] <- t(HX.o.list[[i]])
    HX.u.list[[i]] <- X[H.u.list[[i]], ]
    tHX.u.list[[i]] <- t(HX.u.list[[i]])
  }

  ## initialize spatial covariance
  sigma.squared.eta <- 1 / rgamma(1, alpha.eta, beta.eta)
  phi <- runif(1, phi.lower, phi.upper)
  D.t <- vector('list', length = t)
  Sigma.t <- vector('list', length = t)
  Sigma.t.inv <- vector('list', length = t)
  Sigma.t.star <- vector('list', length = t)
  Sigma.t.inv.star <- vector('list', length = t)
  for(i in 1:t){
    D.t[[i]] <- D[H.list[[i]], H.list[[i]]]
    Sigma.t[[i]] <- I.o[[i]] + sigma.squared.eta * exp( - D.t[[i]] / phi)
    Sigma.t.inv[[i]] <- solve(Sigma.t[[i]])
  }
  
  Sigma.full <- I.full + sigma.squared.eta * exp( - D / phi) 
  Sigma.full.inv <- solve(Sigma.full)
  tX.Sigma.full.inv.X <- tX %*% Sigma.full.inv %*% X
  tX.Sigma.full.inv <- tX %*% Sigma.full.inv

  ## initialize Y.u
  Y.u <- vector('list', length = t)

  projectXontoY <- solve(t(X) %*% Sigma.full.inv %*% X) %*% t(X) %*% Sigma.full.inv
  beta.tilde.gamma <- vector('list', length = t)
  for(i in 1:t){
    if(sum(gamma[[i]]) == 0){
      ##
    } else {
    beta.tilde.gamma[[i]] <- solve(1 / sigma.squared * tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Sigma.t.inv[[i]] %*% HX.o.list[[i]][, gamma[[i]] == 1] + 1 / sigma.squared * Lambda.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Y.list[[i]]
    }
  }
  
  ## initialize sigma.squared
  tmp <- vector(length = t)
  for(i in 1:t){
    if(sum(gamma[[i]]) == 0){
      tmp[i] <- t(Y.list[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
    } else {
      tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
    }
  }
  sigma.squared <- 1 / rgamma(1, (sum(n.o) + sum(unlist(gamma))) / 2, sum(tmp) / 2)
  
  ## initialize variables
  O <- vector('list', length = t)
  rho <- vector('list', length = t)

  ##
  ## setup save variables
  ##
  
  gamma.save <- array(dim = c(p, t, n.mcmc))
  sigma.squared.save <- vector(length = n.mcmc)
  sigma.squared.eta.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc)
  beta.save <- array(dim = c(p, t, n.mcmc))
  rho.save <- array(dim = c(p, t, n.mcmc))
  Y.pred <- array(dim = c(m, t, n.mcmc))
  delta.save <- delta
  phi.accept <- 0
  eta.accept <- 0
  gamma.accept <- 0

  ##
  ## begin mcmc
  ##
  
  for(k in 1:n.mcmc){
#     if(k %% 1000 == 0){
      cat(k, ' ')
#     }
    
    ##
    ## sample Y.u
    ##
    
    for(i in 1:t){
      if(sum(gamma[[i]]) == 0){
        Y.u[[i]] <- Sigma.full[H.u.list[[i]], H.list[[i]]] %*% Sigma.t.inv[[i]] %*% Y.list[[i]]
      } else {
        Y.u[[i]] <- HX.u.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]] + Sigma.full[H.u.list[[i]], H.list[[i]]] %*% Sigma.t.inv[[i]] %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]])
      }
    }
    Y.c <- vector('list', length = t)
    for(i in 1:t){
      Y.c[[i]] <- vector(length = m)
      Y.c[[i]][H.list[[i]]] <- Y.list[[i]]
      Y.c[[i]][H.u.list[[i]]] <- Y.u[[i]]
    }

    beta.hat <- vector('list', length = t)
    projectXontoY <- solve(t(X) %*% Sigma.full.inv %*% X) %*% t(X) %*% Sigma.full.inv
    for(i in 1:t){
      beta.hat[[i]] <- projectXontoY %*% Y.c[[i]]
    }
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      if(sum(gamma[[i]] == 0)){
        tmp[i] <- t(Y.list[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]  
      } else {
        tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% Sigma.t.inv[[i]] %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
      }
      
    }
    sigma.squared <- 1 / rgamma(1, (sum(n.o) + sum(unlist(gamma))) / 2, sum(tmp) / 2)

    ##
    ## sample gammma
    ##

    for(i in 1:t){
      for(j in 1:p){
        if(runif(1) > gamma.tune){
          if(gamma[[i]][j] == 0){
            gamma.star[[i]][j] <- 1
          } else {
            gamma.star[[i]][j] <- 0
          }
        }
      }
    }
    gamma[[1]]
    gamma.star

    mh.gamma.1 <- sum(sapply(1:t, make.gamma.mh, gamma = gamma.star, beta.hat = beta.hat, Sigma.full.inv = Sigma.full.inv, Y.c = Y.c, tX.Sigma.full.inv.X = tX.Sigma.full.inv.X, tX.Sigma.full.inv = tX.Sigma.full.inv))
    mh.gamma.2 <- sum(sapply(1:t, make.gamma.mh, gamma = gamma, beta.hat = beta.hat, Sigma.full.inv = Sigma.full.inv, Y.c = Y.c, tX.Sigma.full.inv.X = tX.Sigma.full.inv.X, tX.Sigma.full.inv = tX.Sigma.full.inv))
    mh.gamma <- exp(mh.gamma.1 - mh.gamma.2)
    
    if(mh.gamma > runif(1)){
      gamma <- gamma.star
      gamma.accept <- 1 / n.mcmc + gamma.accept
    }

# for(i in 1:t){ ## using log scale
#       Psi[[i]] <- 1 / 2 * log(lambda / sigma.squared) - 1 / (2 * sigma.squared) * (beta.hat[[i]]^2 * (lambda - 1000 * delta)) + log(pi.prior) - log(1 - pi.prior)
#       rho[[i]] <- exp(Psi[[i]] - log(1 + exp(Psi[[i]])))
#     }
# 
    for(i in 1:t){
#       gamma[[i]] <- rbinom(p, 1, rho[[i]])
      if(sum(gamma[[i]]) == 0){
        Lambda.gamma[[i]] <- 0
      } else {
      Lambda.gamma[[i]] <- diag(lambda[gamma[[i]] == 1])
      }
    }
    
    ##
    ## sample beta.tilde.gamma
    ##
    
    for(i in 1:t){
      if(sum(gamma[[i]]) == 0){
        beta.tilde.gamma[[i]] <- 0
      } else {
        beta.tilde.gamma[[i]] <- solve(1 / sigma.squared * tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Sigma.t.inv[[i]] %*% HX.o.list[[i]][, gamma[[i]] == 1] + 1 / sigma.squared * Lambda.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Sigma.t.inv[[i]] %*% Y.list[[i]]
      }
    }
    
    ## sample sigma.squared.eta
    sigma.squared.eta.star <- rnorm(1, sigma.squared.eta, sigma.eta.tune)
    if(sigma.squared.eta.star > 0){
      for(i in 1:t){
        Sigma.t.star[[i]] <- I.o[[i]] + sigma.squared.eta.star * exp( - D.t[[i]] / phi)
        Sigma.t.inv.star[[i]] <- solve(Sigma.t.star[[i]])
      }
      mh.eta.1 <- sum(sapply(1:t, make.mh, sigma.squared = sigma.squared, Sigma.t = Sigma.t.star, Sigma.t.inv = Sigma.t.inv.star, gamma = gamma, beta.tilde.gamma = beta.tilde.gamma))
      mh.eta.2 <- sum(sapply(1:t, make.mh, sigma.squared = sigma.squared, Sigma.t = Sigma.t, Sigma.t.inv = Sigma.t.inv, gamma = gamma, beta.tilde.gamma = beta.tilde.gamma))
      mh.eta <- exp(mh.eta.1 - mh.eta.2)
      if(mh.eta > runif(1)){
        sigma.squared.eta <- sigma.squared.eta.star
        Sigma.t <- Sigma.t.star
        Sigma.t.inv <- Sigma.t.inv.star
        eta.accept <- 1 / n.mcmc + eta.accept
      }
    }

    ##
    ## sample phi
    ##

    phi.star <- rnorm(1, phi, phi.tune)
    if(phi.star > phi.lower && phi.star < phi.upper){
      for(i in 1:t){
        Sigma.t.star[[i]] <- I.o[[i]] + sigma.squared.eta * exp( - D.t[[i]] / phi.star)
        Sigma.t.inv.star[[i]] <- solve(Sigma.t.star[[i]])
      }
      mh.phi.1 <- sum(sapply(1:t, make.mh, sigma.squared = sigma.squared, Sigma.t = Sigma.t.star, Sigma.t.inv = Sigma.t.inv.star, gamma = gamma, beta.tilde.gamma = beta.tilde.gamma))
      mh.phi.2 <- sum(sapply(1:t, make.mh, sigma.squared = sigma.squared, Sigma.t = Sigma.t, Sigma.t.inv = Sigma.t.inv, gamma = gamma, beta.tilde.gamma = beta.tilde.gamma))
      mh.phi <- exp(mh.phi.1 - mh.phi.2)
      if(mh.phi > runif(1)){
        phi <- phi.star
        Sigma.t <- Sigma.t.star
        Sigma.t.inv <- Sigma.t.inv.star
        phi.accept <- 1 / n.mcmc + phi.accept
      }
    }

    ##
    ## Sigma.full
    ##
    
    Sigma.full <- I.full * sigma.squared.eta * exp( - D / phi)
    Sigma.full.inv <- solve(Sigma.full)
    tX.Sigma.full.inv.X <- tX %*% Sigma.full.inv %*% X
    tX.Sigma.full.inv <- tX %*% Sigma.full.inv

    ##
    ## log scoring rule
    ##
    
#     log.score <- sum(dnorm(Y.new, mean = cbind(X.new[, 1], X.new[, 2:(p)][, gamma == 1]) %*% beta.tilde.gamma, sd = sqrt(sigma.squared), log = TRUE))
    
    ##
    ## save samples
    ##
    Y.pred[, , k] <- matrix(unlist(Y.c), nrow = m, ncol = t, byrow = FALSE)
    gamma.save[, , k] <- matrix(unlist(gamma), nrow = p, ncol = t, byrow = FALSE)
    sigma.squared.save[k] <- sigma.squared
    sigma.squared.eta.save[k] <- sigma.squared.eta
    phi.save[k] <- phi
    beta.save[, , k] <- matrix(unlist(beta.hat), nrow = p, ncol = t, byrow = FALSE)
#     rho.save[, , k] <- matrix(unlist(rho), nrow = p, ncol = t, byrow = FALSE) 
#     delta.save <- delta
#     log.score.save[k] <- log.score
  }
  list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, rho.save = rho.save, delta.save = delta.save, Y.pred = Y.pred, eta.accept = eta.accept, phi.accept = phi.accept, gamma.accept = gamma.accept, sigma.squared.eta.save = sigma.squared.eta.save, phi.save = phi.save)#, log.score.save = log.score.save)
#   list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, delta.save = delta.save, Y.pred = Y.pred)#, log.score.save = log.score.save)
  
}

