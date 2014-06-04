##
## MCMC sampler for orthogonal data augmentation
##

mcmc.pcaMA <- function(Y.list, X.o, H.list, params, epsilon = 0.001){ #Y.new, X.new, for log scoring rule
  
  ##
  ## functions and subroutines
  ##
  
  ##
  ## initialize fixed values
  ##
  
  n.mcmc <- params[[1]]
  alpha <- params[[2]]
  pi.prior <- params[[3]]  
  lambda <- params[[4]]
  
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
  I.o <- vector('list', length = t)
  I.u <- vector('list', length = t)
  for(i in 1:t){
    I.o[[i]] <- diag(n.o[i])
    I.u[[i]] <- diag(n.u[i])
  }
  ## initialize random values
  
  Psi <- vector('list', length = t)
  gamma <- vector('list', length = t)
  for(i in 1:t){
    gamma[[i]] <- rbinom(p, 1, pi.prior)
  }
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

  ## initialize Y.u
  Y.u <- vector('list', length = t)

  projectXontoY <- solve(t(X) %*% X) %*% t(X)
  beta.tilde.gamma <- vector('list', length = t)
  for(i in 1:t){
    if(sum(gamma[[i]]) == 0){
      ##
    } else {
    beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][gamma[[i]] == 1, ] %*% HX.o.list[[i]][, gamma[[i]] == 1] + Lambda.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Y.list[[i]]
    }
  }
  
  ## initialize sigma.squared
  tmp <- vector(length = t)
  for(i in 1:t){
    tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
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
  beta.save <- array(dim = c(p, t, n.mcmc))
  rho.save <- array(dim = c(p, t, n.mcmc))
  Y.pred <- array(dim = c(m, t, n.mcmc))
  delta.save <- delta

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
        Y.u[[i]] <- rnorm(n.u[i], mean = 0, sd = sigma.squared)
      } else {
        mn = HX.u.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]
        sig.chol = chol(sigma.squared * (I.u[[i]] + HX.u.list[[i]][, gamma[[i]] == 1] %*% solve( 1000 * diag(delta[gamma[[i]] == 1]) + Lambda.gamma[[i]]) %*% tHX.u.list[[i]][gamma[[i]] == 1, ]))
        devs <- rnorm(n.u[i])
        Y.u[[i]] <- sig.chol %*% devs + mn
      }
    }
    Y.c <- vector('list', length = t)
    for(i in 1:t){
      Y.c[[i]] <- vector(length = m)
      Y.c[[i]][H.list[[i]]] <- Y.list[[i]]
      Y.c[[i]][H.u.list[[i]]] <- Y.u[[i]]
    }

    beta.hat <- vector('list', length = t)
    for(i in 1:t){
      beta.hat[[i]] <- projectXontoY %*% Y.c[[i]]
    }
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      if(sum(gamma[[i]] == 0)){
        tmp[i] <- t(Y.list[[i]]) %*% Y.list[[i]] #+ t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
      } else {
        tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Lambda.gamma[[i]] %*% beta.tilde.gamma[[i]]
      }
    }
    sigma.squared <- 1 / rgamma(1, (sum(n.o) + sum(unlist(gamma))) / 2, sum(tmp) / 2)
    
    ##
    ## sample gammma
    ##

    for(i in 1:t){ ## using log scale
      Psi[[i]] <- 1 / 2 * log(lambda / sigma.squared) - 1 / (2 * sigma.squared) * (beta.hat[[i]]^2 * (lambda - 1000 * delta)) + log(pi.prior) - log(1 - pi.prior)
      rho[[i]] <- exp(Psi[[i]] - log(1 + exp(Psi[[i]])))
    }

    for(i in 1:t){
      gamma[[i]] <- rbinom(p, 1, rho[[i]])
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
        ##
      } else {
        beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][gamma[[i]] == 1, ] %*% HX.o.list[[i]][, gamma[[i]] == 1] + Lambda.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Y.list[[i]]
      }
    }
    
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
    beta.save[, , k] <- matrix(unlist(beta.hat), nrow = p, ncol = t, byrow = FALSE)
    rho.save[, , k] <- matrix(unlist(rho), nrow = p, ncol = t, byrow = FALSE) 
#     delta.save <- delta
#     log.score.save[k] <- log.score
  }
  list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, rho.save = rho.save, delta.save = delta.save, Y.pred = Y.pred)#, log.score.save = log.score.save)
#   list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, delta.save = delta.save, Y.pred = Y.pred)#, log.score.save = log.score.save)
  
}

