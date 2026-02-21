#Garch

# Using linear regression
library(quantmod)
library(GIGrvg)
library(fBasics)
library(MASS)

#set.seed(42) # for reproducibility

# Get SPY data for the MCMC
getSymbols("SPY", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
SPY_adj <- Ad(SPY)
R_t <- (SPY_adj / lag(SPY_adj)) - 1
y <- R_t[-1] # Remove NA
y_sd <- sd(as.numeric(y_raw))
n <- length(y)

n_iter <- 300000
burnin <- 50000
thin  <- 100

# Covariates

X <- matrix(1, nrow = length(y), ncol = 1) # Add a column of ones for the intercept for now we only have 1 for intercept
b <- rep(0, length = ncol(X)+1) # Initialize the coefficients vector (0 cause we dont know where it starts)
z <- rep(1, length = length(y)) # Initialize z for the GIG distribution
prior_mean_beta <- rep(0, length = ncol(X)+1) # Prior mean for beta
prior_prec_beta <- diag(1/100, ncol(X)+1) # Prior precision for beta

# Storage Matrices
gamma_store = rep(NA, n_iter) # Initialize the shape parameters
delta_store = rep(NA, n_iter) # Initialize the scale parameter
beta_store <- matrix(NA, nrow = n_iter, ncol = ncol(X)+1) # Initialize the matrix to store the beta coefficients for each iteration, we have ncol(X)+1 because we have the intercept and the slope
z_store <- matrix(NA,nrow=n_iter,ncol=n) # Store z values for each iteration


# Priors

xi <- 0.01 # for gamma φ
x_prior <- 0.01
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
  z <- rgig(n = length(y), lambda = -1, chi = (delta * sqrt(q_val))^2, psi = alpha_param^2) # Sample z from GIG distribution
  
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
  u2 <- n + omega - x_prior
  u3 <- n*z_bar_r + omega/eta
  
  # Use scheme 1 Gamma-GIG
  mu_ig <- rgig(1, lambda = (n-1)/2, chi = phi*u1, psi = phi*u3) # Sample μ from GIG distribution
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
keep_idx   <- seq(burnin + thin, n_iter, by = thin)
beta_post  <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]
z_post_mean <- colMeans(z_store[keep_idx, ])

# Predictive Distribution
n_pred <- 150000
predictions <- numeric(n_pred)

beta_post_median <- apply(beta_post, 2, median)
gamma_post_median <- median(gamma_post)
delta_post_median <- median(delta_post)

for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5, 
                chi = delta_post_median^2, 
                psi = gamma_post_median^2 + beta_post_median[2]^2)
  mu_j  <- beta_post_median[1] + beta_post_median[2] * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}



cat("Tomorrow's SPY return forecast:\n")
cat("Expected:  ", round(mean(predictions)*100, 3), "%\n")
cat("95% CI:    ", round(quantile(predictions, c(0.025,0.975))*100, 3), "%\n")
cat("5% VaR:    ", round(quantile(predictions, 0.05)*100, 3), "%\n")


# plot
par(mfrow = c(1,1), mar = c(5,5,4,2))

hist(y, breaks = 300, probability = TRUE, 
     col = rgb(0.8, 0.9, 1, 0.6), border = "white",
     main = "SPY Daily Returns vs Fitted NIG (Gamma-GIG)",
     xlab = "Daily Return", ylab = "Density",
     ylim = c(0, 75))

lines(density(y, bw = "SJ", adjust = 1.1, n = 1024), 
      col = "black", lwd = 3.2)

x_seq <- seq(min(y)-0.02, max(y)+0.02, length.out = 2000)

lines(x_seq, 
      dnig(x_seq, 
           mu    = beta_post_median[1],
           beta  = beta_post_median[2],
           delta = delta_post_median,
           alpha = sqrt(gamma_post_median^2 + beta_post_median[2]^2)),
      col = "red", lwd = 4)

legend("topright", 
       legend = c("Observed SPY returns", "Fitted NIG"),
       col = c("black", "red"), lwd = c(3.2, 4), bty = "n", cex = 1.15)




# Trace to see convergence
mcmc_samples <- as.mcmc(cbind(beta_post[,1], beta_post[,2], gamma_post, delta_post))
colnames(mcmc_samples) <- c("mu", "beta_skew", "gamma", "delta")
par(mfrow = c(2,2), mar = c(4,4,3,1))
traceplot(mcmc_samples[, "mu"],        main = "Traceplot: μ (location)", ylab = "μ")
traceplot(mcmc_samples[, "beta_skew"], main = "Traceplot: β (skewness)", ylab = "β")
traceplot(mcmc_samples[, "gamma"],     main = "Traceplot: γ", ylab = "γ")
traceplot(mcmc_samples[, "delta"],     main = "Traceplot: δ", ylab = "δ")

# Multiple MCMC
median_mu    <- median(beta_post[, 1])
median_beta  <- median(beta_post[, 2])
median_gamma <- median(gamma_post)
median_delta <- median(delta_post)
median_alpha <- median(sqrt(gamma_post^2 + beta_post[, 2]^2))

cat("=== Posterior Medians (from MCMC) ===\n")
cat("μ     (location)     :", round(median_mu,    6), "\n")
cat("β     (skewness)     :", round(median_beta,  6), "\n")
cat("γ     (shape)        :", round(median_gamma, 5), "\n")
cat("δ     (scale)        :", round(median_delta, 5), "\n")
cat("α     (tail param)   :", round(median_alpha, 5), "\n\n")


