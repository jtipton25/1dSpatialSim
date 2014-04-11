##
## Plot of 1d spatial mcmc output
##

make.output.plot <- function(out){
  n.burn <- floor(n.mcmc / 10)
  #x11()
  layout(matrix(1:9, nrow = 3))
  #
  if(!is.null(out$mu.beta.save)){
    matplot(t(out$mu.beta.save)[(n.burn + 1):n.mcmc, ], type = 'l', main = expression(paste("Trace plot for ", mu[beta])), ylab = expression(mu[beta]), xlab = "MCMC iteration - burn in")
#     abline(h = beta[1], col = 'black')
#     abline(h = beta[2], col = 'red')
  }
  #
  plot(out$sigma.squared.beta.save[(n.burn + 1):n.mcmc], type = 'l', main = expression(paste("Trace plot for ", sigma[beta]^2)), ylab = expression(sigma[beta]^2), xlab = "MCMC iteration - burn in")
  #
  plot(out$sigma.squared.epsilon.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$epsilon.accept, 2)), ylab = expression(sigma[epsilon]^2), xlab = "MCMC iteration - burn in")
  abline(h = s2.e, col = 'red')
  #
  plot(out$sigma.squared.eta.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$eta.accept, 2)), ylab = expression(sigma[eta]^2), xlab = "MCMC iteration - burn in")
  abline(h = s2.s, col = 'red')
  #
  plot(out$phi.save[(n.burn + 1):n.mcmc], type = 'l', main = paste("accept rate", round(out$phi.accept, 2)), ylab = expression(phi), xlab = "MCMC iteration - burn in")
  abline(h = phi, col = 'red')
  #
  matplot(out$fort.raster, type = 'l', ylim = c(min(out$fort.raster) - 2*max(sqrt(out$var.save)), max(out$fort.raster) + 2*max(sqrt(out$var.save))), col = 'black', main = "Fitted Data")
  for(i in 1:dim(out$var.save)[2]){
    matplot(out$fort.raster - 2*sqrt(out$var.save[, i]), type = 'l', add = TRUE, col = adjustcolor('red', alpha = 0.25), lty = 'dashed')
    matplot(out$fort.raster + 2*sqrt(out$var.save[, i]), type = 'l', add = TRUE, col = adjustcolor('red', alpha = 0.25), lty = 'dashed')
  }
  #
  plot.Z.field(Z.list.hist, locs = locs, main = "True Surface")
  #
#   plot.Z.field(Z.list.pca, locs = locs, main = "Pallete of Signals")
  #
  plot.Y.field(field$Y.list, field$H.list, locs = locs)
  #
#   if(!is.null(out$mu.beta.save)){
#     hist(out$mu.beta.save[1, ][(n.burn + 1):n.mcmc])
#     abline(v = beta[1], col = 'red')
#     abline(v = quantile(out$mu.beta.save[1, ], probs = c(0.025, 0.975)), col = 'blue')
#   }
#   #
#   if(!is.null(out$mu.beta.save)){
#     hist(out$mu.beta.save[2, ][(n.burn + 1):n.mcmc])
#     abline(v = beta[2], col = 'red')
#     abline(v = quantile(out$mu.beta.save[2, ], probs = c(0.025, 0.975)), col = 'blue')
#   }
  #
  MSPE <- (out$fort.raster -matrix(unlist(Z.list.hist), nrow = m, byrow = FALSE))^2
  matplot(MSPE, type = 'l', main = paste('MSPE', round(mean(MSPE), 4)), ylab = "MSPE", xlab = "location")
}

