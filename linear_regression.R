# ========================================================
# Robust Bayesian NIG - Gamma-TN Scheme for SPY
# ========================================================

library(quantmod)
library(GIGrvg)
library(fBasics)
library(MASS)

set.seed(42)

# Data
getSymbols("SPY", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
y <- as.numeric(dailyReturn(Ad(SPY))[-1])
n <- length(y)

# Stronger, more stable priors for daily returns
alpha0 <- 8.0      # shape for δ² ~ Gamma
beta0  <- 4.0      # rate for δ²
omega  <- 0.0
theta  <- 3.0      # variance of γ|δ

prior_mean_beta <- c(0, 0)
prior_prec_beta <- diag(1/30, 2)   # slightly stronger prior

# MCMC
n_iter <- 120000
burn   <- 30000
thin   <- 60

beta_store  <- matrix(NA, n_iter, 2)
gamma_store <- numeric(n_iter)
delta_store <- numeric(n_iter)

# Reasonable initial values
mu    <- mean(y)
beta  <- -0.1
gamma <- 40
delta <- sd(y) * 1.5

for (i in 1:n_iter) {
  
  # 1. Latent z
  z <- rgig(n, lambda = -1, chi = (y - mu)^2 + delta^2, psi = gamma^2 + beta^2)
  
  # 2. Update μ and β
  A1 <- cbind(1, z)
  w  <- 1/z
  D_inv <- t(A1 * w) %*% A1 + prior_prec_beta
  D <- solve(D_inv)
  d <- t(A1) %*% (y * w) + prior_prec_beta %*% prior_mean_beta
  b <- mvrnorm(1, D %*% d, D)
  mu   <- b[1]
  beta <- b[2]
  
  # 3. Gamma-TN update
  z_bar   <- mean(z)
  z_bar_r <- mean(1/z)
  
  shape_d <- alpha0 + n/2
  rate_d  <- beta0 + z_bar_r/2 + (omega^2)/(2*theta) - 
    ((omega + n*theta)^2) / (2*theta * (1 + n*z_bar*theta))
  rate_d  <- max(rate_d, 1e-5)
  
  delta2 <- rgamma(1, shape_d, rate_d)
  delta  <- sqrt(max(delta2, 1e-6))
  
  mean_g <- (omega + n*theta) * delta / (1 + n*z_bar*theta)
  sd_g   <- sqrt(theta / (1 + n*z_bar*theta))
  
  gamma <- rnorm(1, mean_g, sd_g)
  if (gamma <= 0) gamma <- abs(mean_g) * 0.5 + 0.05
  
  # Store
  beta_store[i,]  <- c(mu, beta)
  gamma_store[i]  <- gamma
  delta_store[i]  <- delta
}

# ====================== POSTERIOR ======================
keep <- seq(burn + 1, n_iter, by = thin)
beta_post  <- beta_store[keep, ]
gamma_post <- gamma_store[keep]
delta_post <- delta_store[keep]

beta_post_median <- apply(beta_post, 2, median)
gamma_post_median <- median(gamma_post)
delta_post_median <- median(delta_post)

cat("Posterior Medians:\n")
cat("μ     :", round(beta_post_median[1], 6), "\n")
cat("β     :", round(beta_post_median[2], 6), "\n")
cat("γ     :", round(gamma_post_median, 4), "\n")
cat("δ     :", round(delta_post_median, 5), "\n")
cat("α     :", round(sqrt(gamma_post_median^2 + beta_post_median[2]^2), 4), "\n\n")

# ====================== PLOT (Normal scale - should look good now) ======================
par(mfrow = c(1,1), mar = c(5,5,4,2))

hist(y, breaks = 250, probability = TRUE, col = rgb(0.8,0.9,1,0.7), border = "white",
     main = "SPY Daily Returns vs Fitted NIG (Gamma-TN)",
     xlab = "Daily Return", ylab = "Density")

x_seq <- seq(min(y)-0.02, max(y)+0.02, length.out = 2000)
lines(x_seq, dnig(x_seq, 
                  mu = beta_post_median[1],
                  beta = beta_post_median[2],
                  delta = delta_post_median,
                  alpha = sqrt(gamma_post_median^2 + beta_post_median[2]^2)),
      col = "red", lwd = 3.5)

lines(density(y, bw = "SJ", adjust = 1.1), col = "black", lwd = 3)

legend("topright", legend = c("Observed SPY returns", "Fitted NIG (Gamma-TN)"),
       col = c("black", "red"), lwd = 3, bty = "n")