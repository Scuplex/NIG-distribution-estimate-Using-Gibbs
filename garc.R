#Garch

# Using linear regression
library(quantmod)
library(GIGrvg)
library(fBasics)
library(MASS)
library(moments)
library(tseries)
set.seed(42) # for reproducibility

# Get SPY data for the MCMC
getSymbols("SPY", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
SPY_adj <- Ad(SPY)
R_t <- (SPY_adj / lag(SPY_adj)) - 1
y <- R_t[-1] # Remove NA
n <- length(y)
n_iter <- 100000
burnin <- 20000

# Covariates
getSymbols("DGS10", src = "FRED", from = "2014-12-01", to = "2025-12-31")
yield_change <- diff(DGS10)
yield_change <- na.omit(yield_change)
common_dates <- intersect(index(yield_change), index(R_t[-1]))
yield_aligned <- as.vector(yield_change[common_dates])
y <- as.vector(R_t[-1][common_dates])
n <- length(y)
ADF_result <- adf.test(y)
cat("ADF Test p-value:", ADF_result$p.value, "\n")


# Matrix X and Y

#X <- matrix(1, nrow = length(y), ncol = 1) # Add a column of ones for the intercept for now we only have 1 for intercept
X <- cbind(1, yield_aligned)
beta_store <- matrix(NA, nrow = n_iter, ncol = ncol(X)+1) # Initialize the matrix to store the beta coefficients for each iteration, we have ncol(X)+1 because we have the intercept and the slope
b <- rep(0, length = ncol(X)+1) # Initialize the coefficients vector (0 cause we dont know where it starts)
prior_mean_beta <- rep(0, length = ncol(X)+1) # Prior mean for beta
prior_prec_beta <- diag(1/100, ncol(X)+1) # Prior precision for beta
gamma_store = rep(NA, n_iter) # Initialize the shape parameters
delta_store = rep(NA, n_iter) # Initialize the scale parameter
z <- rep(1, length = length(y)) # Initialize z for the GIG distribution
z_store <- matrix(NA,nrow=n_iter,ncol=n) # Store z values for each iteration

xi <- 0.01 # for gamma φ
x_prior <- 0.01 # for gamma φ

phi <- 1 # for IG distribution
eta <- var(y) # prior mean for μ in the IG distribution
omega <- 0.001 # Prior how strongly we believe in the prior mean for μ in the IG distribution

gamma <- 1 
delta <- 1 
mu_ig <- delta/gamma # prior mean for the μ parameter in the IG distribution


for (i in 1:n_iter) {
  
  # Z draws
  residuals <- y - X %*% b[1:ncol(X)] # Calculate residuals y-μ1
  q_val <- 1 + (residuals/delta)^2
  alpha_param <- sqrt(gamma^2 + b[ncol(X)+1]^2) # alpha parameter
  z <- rgig(n = length(y), lambda = -1, chi = delta * sqrt(q_val), psi = alpha_param) # Sample z from GIG distribution
  
  # d,D matrices and coefficient generation
  A1    <- cbind(X, z)                    # n x (k+1) design matrix
  w <- as.vector(1/z)            # n x n precision vector diag(1/z1,...,1/zn)
  A1w <- A1 * w                  # n x (k+1) weighted design matrix for less computation
  C2inv <- prior_prec_beta                # (k+1) x (k+1) prior precision matrix
  
  D_inv <- t(A1w) %*% A1 + C2inv # (k+1) x (k+1) posterior precision matrix
  D     <- solve(D_inv) # (k+1) x (k+1) posterior covariance matrix
  d <- t(A1) %*% (y*w) + C2inv %*% prior_mean_beta # posterior mean
  
  b     <- mvrnorm(1, mu = D %*% d, Sigma = D) # draw 1 sample of beta from posterior MVN(D*d, D)
  
  # Find statistic values to compute gamma delta
  z_bar <- mean(z) 
  z_bar_r <- mean(1/z) 
  nu <- n + 2*xi 
  u1 <- n*z_bar + omega*eta 
  u2 <- n + omega - xi
  u3 <- n*z_bar_r + omega/eta
  
  # Use scheme 1 Gamma-GIG
  mu_ig <- rgig(1, lambda =(n-1)/2, chi = sqrt(phi*u1), psi=sqrt(phi*u3)) # Sample μ from GIG distribution
  phi <- rgamma(1,shape=(nu+1)/2,rate=u1/(2*mu_ig)-u2+u3*mu_ig/2) # Sample φ from gamma distribution
  
  #Convert to Gamma,Delta
  gamma <- sqrt(phi / mu_ig)
  delta <- sqrt(mu_ig * phi)
  
  beta_store[i,] <- b # Store the sampled beta coefficients
  gamma_store[i] <- gamma # Store the sampled gamma values
  delta_store[i] <- delta # Store the sampled delta values
  z_store[i,] <- z # Store the sampled z values
  
}

# Discard burn-in samples
thin       <- 10
keep_idx   <- seq(burnin + thin, n_iter, by = thin)
beta_post  <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]
z_post_mean <- colMeans(z_store[keep_idx, ])

# Predictive Distribution
#Assum Yield_tomorrow = 0 ( we have to do a regression model for it )

yield_tomorrow <- 0  # assume no change
n_sim <- nrow(beta_post)
scale_factor <- var(y) / mean(z_post_mean)

predictions <- rep(NA, n_sim)
for(j in 1:n_sim){
  z_new <- rgig(1, lambda=-1, chi=delta_post[j]^2, psi=gamma_post[j]^2 + beta_post[j,3]^2)
  z_scaled <- z_new * scale_factor
  mu_j <- beta_post[j,1] + beta_post[j,2]*0 + beta_post[j,3]*z_scaled
  predictions[j] <- rnorm(1, mu_j, sqrt(z_scaled))
}

cat("Tomorrow's SPY return forecast:\n")
cat("Expected:  ", round(mean(predictions)*100, 3), "%\n")
cat("95% CI:    ", round(quantile(predictions, c(0.025,0.975))*100, 3), "%\n")
cat("5% VaR:    ", round(quantile(predictions, 0.05)*100, 3), "%\n")

# # Trace Plots

# par(mfrow=c(2,2))
# plot(beta_post[,1], type='l', main='Trace Plot for Intercept', xlab='Iteration', ylab='Value')
# plot(beta_post[,2], type='l', main='Trace Plot for Slope', xlab='Iteration', ylab='Value')
# plot(gamma_post, type='l', main='Trace Plot for Gamma', xlab='Iteration', ylab='Value')
# plot(delta_post, type='l', main='Trace Plot for Delta', xlab='Iteration', ylab='Value')
# 
# # Plot zi over time
# par(mfrow=c(1,1))
# plot(z_post_mean, type='l', main='Mean Z values over time', xlab='Time', ylab='Mean Z value')
# 
# 
# # Posterior Density Plots with statistics
# par(mfrow=c(2,2))
# 
# # Intercept
# plot(density(beta_post[,1]), main='Posterior Density for Intercept', xlab='Value', ylab='Density')
# abline(v=mean(beta_post[,1]), col='red', lwd=2)
# abline(v=quantile(beta_post[,1], c(0.025, 0.975)), col='blue', lwd=2, lty=2)
# 
# # Slope
# plot(density(beta_post[,2]), main='Posterior Density for Slope', xlab='Value', ylab='Density')
# abline(v=mean(beta_post[,2]), col='red', lwd=2)
# abline(v=quantile(beta_post[,2], c(0.025, 0.975)), col='blue', lwd=2, lty=2)
# 
# # Gamma
# plot(density(gamma_post), main='Posterior Density for Gamma', xlab='Value', ylab='Density')
# abline(v=mean(gamma_post), col='red', lwd=2)
# abline(v=quantile(gamma_post, c(0.025, 0.975)), col='blue', lwd=2, lty=2)
# 
# # Delta
# plot(density(delta_post), main='Posterior Density for Delta', xlab='Value', ylab='Density')
# abline(v=mean(delta_post), col='red', lwd=2)
# abline(v=quantile(delta_post, c(0.025, 0.975)), col='blue', lwd=2, lty=2)
# 
# # Summary
# 
# # Summary for Regression Coefficients (Intercept and Beta/Skewness)
# # We calculate Mean, and the 95% Credible Interval (2.5% and 97.5%)
# beta_summary <- data.frame(
#   Mean = colMeans(beta_post),
#   Lower_95 = apply(beta_post, 2, quantile, probs = 0.025),
#   Upper_95 = apply(beta_post, 2, quantile, probs = 0.975)
# )
# rownames(beta_summary) <- c("Intercept", "Yield Change", "Skewness (Beta)")
# 
# # Summary for NIG Scale and Shape
# nig_summary <- data.frame(
#   Parameter = c("Gamma (Shape)", "Delta (Scale)"),
#   Mean = c(mean(gamma_post), mean(delta_post)),
#   Lower_95 = c(quantile(gamma_post, 0.025), quantile(delta_post, 0.025)),
#   Upper_95 = c(quantile(gamma_post, 0.975), quantile(delta_post, 0.975))
# )
# 
# print(beta_summary)
# print(nig_summary)





