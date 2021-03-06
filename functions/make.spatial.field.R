##
## Libraries and Subroutines
##

## make true spatial field
## change the nugget effect
# make.Z.list <- function(reps, mu, Sig.s, m){
#   mu + t(chol(Sig.s)) %*% rnorm(m)
# }

make.Z.list <- function(reps, mu, Sig, m){
  mu + t(chol(Sig)) %*% rnorm(m)
}

## make sampling matrix H
make.H.list <- function(reps, samp, m){
  (1:m)[samp[[reps]]]
}

## change the nugget effect
## make sample data Y
# make.Y.list <- function(reps, Z.list, H.list, s2.e){
#   Z.list[[reps]][H.list[[reps]]] + rnorm(length(H.list[[reps]]), s2.e)
# }

make.Y.list <- function(reps, Z.list, H.list){
  Z.list[[reps]][H.list[[reps]]]
}

plot.Z.field <- function(Z.list, locs, main = "Observed Data", ylab = "Y", xlab = "X"){
  reps <- length(Z.list)
  min.Z <- min(unlist(lapply(Z.list, min)))
  max.Z <- max(unlist(lapply(Z.list, max)))
  plot(Z.list[[1]] ~ locs, type = 'l', ylim = c(min.Z, max.Z), main = main, ylab = ylab, xlab = xlab, col = 1)
  for(t in 2:reps){
    lines(Z.list[[t]] ~ locs, type = 'l', col = t)
  }
}

plot.Y.field <- function(Y.list, H.list, locs, main = "Observed Data", ylab = "Y", xlab = "X"){
  reps <- length(Y.list)
  min.Y <- min(unlist(lapply(Y.list, min)))
  max.Y <- max(unlist(lapply(Y.list, max)))
  idx <- order(locs[H.list[[1]]])
  plot(Y.list[[1]][idx] ~ locs[H.list[[1]]][idx], type = 'l', ylim = c(min.Y, max.Y), main = main, ylab = ylab, xlab = xlab)
  for(t in 2:reps){
    idx <- order(locs[H.list[[t]]])
    lines(Y.list[[t]][idx] ~ locs[H.list[[t]]][idx], type = 'l', col = t)
  }
}

####
####  Simulate 1-D spatial random fields with trend
####

make.spatial.field <- function(reps, X, beta, locs, param = c(s2.s, phi), method = 'exponential', s2.e, samp.size){
  mu <- X %*% beta # mean function
  ## Exponential Spatial Decay Function s2.s * exp( - D / phi)
  if(method == 'exponential'){
    s2.s <- param[1]
    phi <- param[2]
    D <- as.matrix(dist(locs)) # distance matrix
    Sig.s <- s2.s * exp( - D / phi) # spatial covariance matrix
    Sig <- Sig.s + s2.e * diag(dim(D)[2])
    Sig.s.inv <- solve(Sig.s) 
  }

  ## Simulate Random Field with nugget
#   Z.list <- lapply(1:reps, make.Z.list, mu = mu, Sig.s = Sig.s, m = m)
    Z.list <- lapply(1:reps, make.Z.list, mu = mu, Sig = Sig, m = m)
           
  ##  Subsample Fields    
  #samp <- rep(list(sample(1:m, samp.size)), reps)
  tmp.vec = sin(0.025 * (1:m))
  p.vec = tmp.vec - min(tmp.vec)
  samp <- vector('list', length = reps)
  for(i in 1:reps){
      samp[[i]] <- sample(1:m, sample(samp.size, 1), prob = p.vec)
  }
  H.list <- lapply(1:reps, make.H.list, samp = samp, m = m)
#   Y.list <- lapply(1:reps, make.Y.list, Z.list = Z.list, H.list = H.list, s2.e = s2.e)
    Y.list <- lapply(1:reps, make.Y.list, Z.list = Z.list, H.list = H.list)
  ## write output
  list(Z.list = Z.list, Y.list = Y.list, H.list = H.list)
}
