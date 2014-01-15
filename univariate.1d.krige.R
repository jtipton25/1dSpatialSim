##
## Creates data for use in mcmc.spatial kriging code
##

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
make.field <- function(s, mu, Sig.s, m){
  mu + t(chol(Sig.s)) %*% rnorm(m)
}

## make sampling matrix H.list
make.H.list <- function(t, samp, m){
  (1:m)[samp[[t]]]
}

## make sample data Y.list
make.Y.list <- function(t, Z.list, H.list, s2.e = 1){
  Z.list[[t]][H.list[[t]]] + rnorm(length(H.list[[t]]), s2.e)
}

## make sample spatial covaraites X.list
make.X.list <- function(s, X, H.list){
  X[H.list[[s]], ]
}

plot.field <- function(Y.list, H.list, locs){
  t <- length(Y.list)
  min.Y <- min(unlist(lapply(Y.list, min)))
  max.Y <- max(unlist(lapply(Y.list, max)))
  idx <- order(locs[H.list[[1]]])
  plot(Y.list[[1]][idx] ~ locs[H.list[[1]]][idx], type = 'l', ylim = c(min.Y, max.Y), main = "Observed Data")
  if(t > 1){
    for(i in 2:t){
      idx <- order(locs[H.list[[i]]])
      lines(Y.list[[i]][idx] ~ locs[H.list[[i]]][idx], type = 'l', col = i)
    }
  }
}

####
####  Simulate 1-D spatial random fields with trend
####

m <- 1000 # number of spatial locations
locs <- seq(0, 1, , m) # spatial coordinate
t <- 1 # number of spatial fields

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

Z.list <- lapply(1:t, make.field, mu = mu, Sig.s = Sig.s, m = m)

##
## Image full data
##

plot.field(Z.list, H.list = list(rep(1:length(Z.list[[1]]))), locs = locs)

####
####  Subsample Fields and Create Data Matrix   
####

s2.e <- 0.1 # sampling error
# samp.size <- sample(140:200, t, replace = TRUE) # sample size varies
samp.size <- sample(40:200, 1, replace = TRUE) # sample size consistent <- needed for mcmc.spatial

samp <- vector('list', length = t)
for(i in 1:t){
  #samp[[i]] <- sample(1:m, samp.size[[i]])
  samp[[i]] <- sample(1:m, samp.size)
}

H.list <- lapply(1:t, make.H.list, samp = samp, m = m)
Y.list <- lapply(1:t, make.Y.list, Z.list = Z.list, H.list = H.list, s2.e = s2.e)
X.list <- lapply(1:t, make.X.list, X, H.list)

##
## Plot Sample Data
##

plot.field(Y.list = Y.list, H.list = H.list, locs = locs)
#save.image(file = "1dKrige.RData")
