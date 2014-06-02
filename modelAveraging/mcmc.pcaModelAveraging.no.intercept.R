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
  
  gamma <- vector('list', length = t)
  for(i in 1:t){
    gamma[[i]] <- rbinom(p, 1, pi.prior)
  }
  Delta.gamma <- vector('list', length = t)
  for(i in 1:t){
    if(sum(gamma[[i]]) == 0){
      Delta.gamma[[i]] <- 0
    } else {
    Delta.gamma[[i]] <- diag(lambda[which(gamma[[i]] == 1)])  
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
    beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][gamma[[i]] == 1, ] %*% HX.o.list[[i]][, gamma[[i]] == 1] + Delta.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Y.list[[i]]
  }
  
  ## initialize sigma.squared
  tmp <- vector(length = t)
  for(i in 1:t){
    tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Delta.gamma[[i]] %*% beta.tilde.gamma[[i]]
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
#   delta.save <- delta
#   log.score.save <- vector(length = n.mcmc)  

  ##
  ## begin mcmc
  ##
  
  for(k in 1:n.mcmc){
#     if(k %% 10000 == 0){
      cat(k, ' ')
#     }
    
    ##
    ## sample Y.u
    ##
    
    for(i in 1:t){
      if(sum(gamma[[i]] == 0)){
        Y.u[[i]] <- rnorm(n.u[i], mean = 0, sd = sigma.squared)
      } else {
        Y.u[[i]] <- rnorm(n.u[i], mean = HX.u.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]], sd = sigma.squared)
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
        tmp[i] <- t(Y.list[[i]]) %*% Y.list[[i]] #+ t(beta.tilde.gamma[[i]]) %*% Delta.gamma[[i]] %*% beta.tilde.gamma[[i]]
      } else {
        tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, gamma[[i]] == 1] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Delta.gamma[[i]] %*% beta.tilde.gamma[[i]]
      }
    }
    sigma.squared <- 1 / rgamma(1, (sum(n.o) + sum(unlist(gamma))) / 2, sum(tmp) / 2)
    
    ##
    ## sample gammma
    ##

    for(i in 1:t){
      O[[i]] <- log(pi.prior / (1 - pi.prior) * (lambda / (delta + lambda))^(1 / 2)) + 1 / 2 * delta / (delta + lambda * beta.hat[[i]]^2 / sigma.squared * delta)
      rho[[i]] <- exp(O[[i]]) / (1 + exp(O[[i]]))
      for(j in 1:p){
        if(is.na(rho[[i]][j])){
          rho[[i]][j] <- 1
        }
      }
    }
    
    for(i in 1:t){
      gamma[[i]] <- rbinom(p, 1, rho[[i]])
      if(sum(gamma[[i]]) == 0){
        Delta.gamma[[i]] <- 0
      } else {
      Delta.gamma[[i]] <- diag(lambda[gamma[[i]] == 1])
      }
    }
    
    ##
    ## sample beta.tilde.gamma
    ##
    
    for(i in 1:t){
      if(sum(gamma[[i]] == 0)){
        ##
      } else {
        beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][gamma[[i]] == 1, ] %*% HX.o.list[[i]][, gamma[[i]] == 1] + Delta.gamma[[i]]) %*% tHX.o.list[[i]][gamma[[i]] == 1, ] %*% Y.list[[i]]
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
    delta.save <- delta
#     log.score.save[k] <- log.score
  }
  list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, rho.save = rho.save, delta.save = delta.save, Y.pred = Y.pred)#, log.score.save = log.score.save)
  
}

