rm(list = ls())
set.seed(201)

##
## Libraries and Subroutines
##
library('geoR')

## make true spatial field
make.field <- function(reps, mu, Sig.s, m){
	mu + t(chol(Sig.s)) %*% rnorm(m)
}

## make sampling matrix H
make.H.list <- function(reps, samp, m){
  (1:m)[samp[[reps]]]
}

## make sample data Y
make.Y.list <- function(reps, Z.list, H.list, s2.e = 1){
  Z.list[[reps]][H.list[[reps]]] + rnorm(length(H.list[[reps]]), s2.e)
}

make.krige.fit <- function(reps, Y.list, H.list, locs){
  Y.geo <- as.geodata(cbind(Y.list[[reps]], locs[H.list[[reps]]], 
    rep(0, length(H.list[[reps]]))), coords.col = 2:3, data.col = 1)
  
  Y.geo.var <- variog(Y.geo, ini.cov.pairs = c(1, 2), trend = '1st', 
    max.dist = 0.75, messages = FALSE)
  
  Y.geo.fit <- variofit(Y.geo.var, messages = FALSE, cov.model = 'exponential', 
    weights = 'cressie', fix.nugget = TRUE, nugget = 0, fix.kappa = TRUE, kappa = 0) 
  
  return(Y.geo.fit)
}

make.pred <- function(reps, Y.list, H.list, mu, Sig){
	pred <- vector(length = m)
	pred[ - H.list[[reps]]] <- mu[ - H.list[[reps]]] + Sig[[reps]][ - H.list[[reps]], H.list[[reps]]] %*% 
		solve(Sig[[reps]][H.list[[reps]], H.list[[reps]]]) %*% (Y.list[[reps]] - mu[H.list[[reps]]]) 
	pred[H.list[[reps]]] <- Y.list[[reps]]
	return(pred)
}

plot.field <- function(Y.list, H.list, locs, main = "Observed Data", ylab = "Y", xlab = "X"){
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

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
reps <- 100 # number of spatial fields

X <- cbind(rep(1, m), locs)
beta <- c(0, 2) # beta
mu <- X %*% beta # mean function

##
## Exponentail Spatial Decay Function s2.s * exp( - D / phi)
##

s2.s <- 1 # spatial variance parameter
phi <- 15 # spatial decay for exponential
D <- as.matrix(dist(locs)) # distance matrix
Sig.s <- s2.s * exp( - D / phi) # spatial covariance matrix
Sig.s.inv <- solve(Sig.s) 

##
## Simulate Random Field without nugget
##

Z.list <- lapply(1:reps, make.field, mu = mu, Sig.s = Sig.s, m = m)

##
## Image full data
##

H.list.full <-rep(list(1:m), reps)
plot.field(Z.list, H.list = H.list.full, locs = locs, main = "Full Data")
           
####
####  Subsample Fields and Create Data Matrix   
####

s2.e <- 0.01 # sampling error
#samp.size <- sample(140:200, reps, replace = TRUE) # sample size
samp.size <- sample(100:200, 1)

samp <- rep(list(sample(1:m, samp.size)), reps)
#samp <- vector('list', length = reps)

#for(t in 1:reps){
#  samp[[t]] <- sample(1:m, samp.size[[t]])
#}

H.list <- lapply(1:reps, make.H.list, samp = samp, m = m)
Y.list <- lapply(1:reps, make.Y.list, Z.list = Z.list, H.list = H.list, s2.e = s2.e)

plot.field(Y.list = Y.list, H.list = H.list, locs = locs, main = "Sampled Data")











#save.image(file = "1dKrige.RData")

##
## Krige the full data
##

## Assume stationarity and anisotropy
Z.geo <- list(length = reps)
Z.geo.var <- list(length = reps)
Z.geo.fit <- list(length = reps)

for(t in 1:reps){
  Z.geo[[t]] <- as.geodata(cbind(Z.list[[t]], locs, rep(0, m)), coords.col = 2:3, data.col = 1)
  Z.geo.var[[t]] <- variog(Z.geo[[t]], ini.cov.pairs = c(1, 2), trend = '1st', 
                           max.dist = 0.75, messages = FALSE)
  Z.geo.fit[[t]] <- variofit(Z.geo.var[[t]], messages = FALSE, cov.model = 'exponential', 
                             weights = 'cressie', fix.nugget = TRUE, nugget = 0, fix.kappa = TRUE, kappa = 0) 
}

## Plot Semivariogram
for(t in 1:reps){
  plot(Z.geo.var[[t]])
  lines(Z.geo.fit[[t]])
}

s2.fit <- vector(length = reps)
phi.fit <- vector(length = reps)
Sig.s.fit <- vector('list', length = reps)

for(t in 1:reps){
  s2.fit[t] <- summary(Z.geo.fit[[t]])$spatial.component[1]
  phi.fit[t] <- summary(Z.geo.fit[[t]])$spatial.component[2]
  Sig.s.fit[[t]] <- s2.fit[t] * exp( - D / phi.fit[t])
}


##
## Krige the sample data
##

##  coords.col = 2:3, data.col = 1)

##
## Predict spatial field using Kriging with known covariance
##

Sig <- vector('list', length = reps)

for(t in 1:reps){
  Sig[[t]] <- Sig.s
}

Z.pred <- lapply(1:reps, make.pred, Y.list = Y.list, H.list = H.list, mu = mu, Sig = Sig)
Z.pred <- matrix(unlist(Z.pred), nrow = m)
matplot(Z.pred, type = 'l')

##
## Krige using covariance matrix estimated from sample data
##

## Has problems due to struggles estimating phi

Y.geo.fit <- lapply(1:reps, make.krige.fit, Y.list = Y.list, H.list = H.list, locs = locs)

s2.fit <- vector(length = reps)
phi.fit <- vector(length = reps)
Sig.s.fit <- vector('list', length = reps)

for(t in 1:reps){
	s2.fit[t] <- summary(Y.geo.fit[[t]])$spatial.component[1]
	phi.fit[t] <- summary(Y.geo.fit[[t]])$spatial.component[2]
	Sig.s.fit[[t]] <- s2.fit[t] * exp( - D / phi.fit[t])
}

## phi is currently underestimated

## Predict spatial field Using Kriging with estimated covariance
Z.pred.knownS <- lapply(1:reps, make.pred, Y.list)
mu.hat <- lapply(Y.list, mean)

Z.pred <- lapply(1:reps, make.pred, Y.list = Y.list, H.list = H.list, mu = mu, Sig = Sig.s.fit)
  for(t in 1:reps){
Z.pred <- vector(length = m)
Z.pred[ - H.list] <- mu[ - H.list] + Sig.s[ - H.list, H.list] %*% 
  solve(Sig.s[H.list, H.list]) %*% (Y.list - mu[H.list]) 
Z.pred[H.list] <- Y.list
}

## predice spatial field using estimated mean and covariance

##
## Image full data and predicted surface
##

layout(matrix(1:2, 2))
matplot(Z.pred.knownS, type = 'l', lty = 1, main = "Predicted 1-D Spatial Field")
matplot(Z, type = "l", lty = 1, main = "1-D Spatial Field")

