rm(list = ls())

##
## Libraries and Subroutines
##

library(geoR)
#library(gstat)
library(mvtnorm)
library(gplots)
set.seed(101)

## make true spatial field
make.field <- function(t, mu, Sig.s, m){
	mu + t(chol(Sig.s)) %*% rnorm(m)
}

## make sampling matrix H.list
make.H.list <- function(t, samp, m){
  (1:m)[samp[[t]]]
}

## make sample data Y.list
make.Y.list <- function(t, Z, H.list, s2.e = 1){
  Z[, t][H.list[[t]]] + rnorm(length(H.list[[t]]), s2.e)
}

## make sample spatial covaraites X.list
make.X.list <- function(s, X, H.list){
  X[H.list[[s]], ]
}


make.krige.fit <- function(t, Y.list, H.list, locs){
  Y.geo <- as.geodata(cbind(Y.list[[t]], locs[H.list[[t]]], 
    rep(0, length(H.list[[t]]))), coords.col = 2:3, data.col = 1)
  
  Y.geo.var <- variog(Y.geo, ini.cov.pairs = c(1, 2), trend = '1st', 
    max.dist = 0.75, messages = FALSE)
  
  Y.geo.fit <- variofit(Y.geo.var, messages = FALSE, cov.model = 'exponential', 
    weights = 'cressie', fix.nugget = TRUE, nugget = 0, fix.kappa = TRUE, kappa = 0) 
  
  return(Y.geo.fit)
}

make.pred <- function(t, Y.list, H.list, mu, Sig){
	pred <- vector(length = m)
	pred[ - H.list[[t]]] <- mu[ - H.list[[t]]] + Sig[[t]][ - H.list[[t]], H.list[[t]]] %*% 
		solve(Sig[[t]][H.list[[t]], H.list[[t]]]) %*% (Y.list[[t]] - mu[H.list[[t]]]) 
	pred[H.list[[t]]] <- Y.list[[t]]
	return(pred)
}

plot.field <- function(Y.list, H.list, locs){
  t <- length(Y.list)
  min.Y <- min(unlist(lapply(Y.list, min)))
  max.Y <- max(unlist(lapply(Y.list, max)))
  idx <- order(locs[H.list[[1]]])
  plot(Y.list[[1]][idx] ~ locs[H.list[[1]]][idx], type = 'l', ylim = c(min.Y, max.Y), main = "Observed Data")
  for(i in 2:t){
    idx <- order(locs[H.list[[i]]])
    lines(Y.list[[i]][idx] ~ locs[H.list[[i]]][idx], type = 'l', col = i)
  }
}

####
####  Simulate 1-D spatial random fields with trend
####

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
t <- 100 # number of spatial fields

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

Z <- sapply(1:t, make.field, mu = mu, Sig.s = Sig.s, m = m)

##
## Image full data
##

matplot(Z, type = "l", lty = 1, main = "1-D Spatial Field")

##
## Krige the full data
##

## Assume stationarity and anisotropy
Z.geo <- list(length = t)
Z.geo.var <- list(length = t)
Z.geo.fit <- list(length = t)

for(i in 1:t){
  Z.geo[[i]] <- as.geodata(cbind(Z[, i], locs, rep(0, m)), coords.col = 2:3, data.col = 1)
  Z.geo.var[[i]] <- variog(Z.geo[[i]], ini.cov.pairs = c(1, 2), trend = '1st', 
    max.dist = 0.75, messages = FALSE)
  Z.geo.fit[[i]] <- variofit(Z.geo.var[[i]], messages = FALSE, cov.model = 'exponential', 
    weights = 'cressie', fix.nugget = TRUE, nugget = 0, fix.kappa = TRUE, kappa = 0) 
}

## Plot Semivariogram
for(i in 1:t){
  plot(Z.geo.var[[i]])
  lines(Z.geo.fit[[i]])
}

s2.fit <- vector(length = t)
phi.fit <- vector(length = t)
Sig.s.fit <- vector('list', length = t)
	
for(i in 1:t){
  s2.fit[i] <- summary(Z.geo.fit[[i]])$spatial.component[1]
  phi.fit[i] <- summary(Z.geo.fit[[i]])$spatial.component[2]
  Sig.s.fit[[i]] <- s2.fit[i] * exp( - D / phi.fit[i])
}

####
####  Subsample Fields and Create Data Matrix   
####

s2.e <- 0.01 # sampling error
# samp.size <- sample(140:200, t, replace = TRUE) # sample size varies
samp.size <- sample(40:200, 1, replace = TRUE) # sample size consistent <- needed for mcmc.spatial

samp <- vector('list', length = t)
for(i in 1:t){
  #samp[[i]] <- sample(1:m, samp.size[[i]])
  samp[[i]] <- sample(1:m, samp.size)
}

H.list <- lapply(1:t, make.H.list, samp = samp, m = m)
Y.list <- lapply(1:t, make.Y.list, Z = Z, H.list = H.list, s2.e = s2.e)
X.list <- lapply(1:t, make.X.list, X, H.list)

##
## Plot Sample Data
##

plot.field(Y.list = Y.list, H.list = H.list, locs = locs)

#save.image(file = "1dKrige.RData")

##
## Krige the data
##

##  coords.col = 2:3, data.col = 1)

##
## Predict spatial field using Kriging with known covariance
##

Sig <- vector('list', length = t)

for(i in 1:t){
  Sig[[i]] <- Sig.s
}

Z.pred <- lapply(1:t, make.pred, Y.list = Y.list, H.list = H.list, mu = mu, Sig = Sig)
Z.pred <- matrix(unlist(Z.pred), nrow = m)
matplot(Z.pred, type = 'l')

##
## Krige using covariance matrix estimated from sample data
##

## Has problems due to struggles estimating phi

Y.geo.fit <- lapply(1:t, make.krige.fit, Y.list = Y.list, H.list = H.list, locs = locs)

s2.fit <- vector(length = t)
phi.fit <- vector(length = t)
Sig.s.fit <- vector('list', length = t)

for(i in 1:t){
	s2.fit[i] <- summary(Y.geo.fit[[i]])$spatial.component[1]
	phi.fit[i] <- summary(Y.geo.fit[[i]])$spatial.component[2]
	Sig.s.fit[[i]] <- s2.fit[i] * exp( - D / phi.fit[t])
}

## phi is currently underestimated

## Predict spatial field Using Kriging with estimated covariance
Z.pred.knownS <- lapply(1:t, make.pred, Y.list)
mu.hat <- lapply(Y.list, mean)

Z.pred <- lapply(1:t, make.pred, Y.list = Y.list, H.list = H.list, mu = mu, Sig = Sig.s.fit)
  for(i in 1:t){
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

