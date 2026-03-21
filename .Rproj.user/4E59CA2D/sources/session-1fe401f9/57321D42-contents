# Bayesian estimation of heteroscedastic regression

library(quantmod)
library(GIGrvg)
library(MASS)
library(parallel)
library(fBasics)
set.seed(42) # for reproducibility
mem.maxVSize(64000)

# Get SPY data for the MCMC

ticker <- "^SP500TR"
start_date <- "2015-12-01"
end_date <- "2025-12-31"

data_raw <- getSymbols(ticker, from = start_date, to = end_date, periodicity = "daily", auto.assign = FALSE)
y <- 100 * as.numeric(dailyReturn(Ad(data_raw), type = "arithmetic"))[-1]  # ← percent returns
n <- length(y)

N_iteration <- 500000
burnin <- 150000
thin  <- 200
n_pred <- 50000


# Initialize Variables/Matrices
X <- matrix(1, nrow = n, ncol = 1) # Add a column of ones for the intercept for now we only have 1 for intercept
b <- rep(0, ncol(X)+1) # Initialize the coefficients vector (0 cause we don't know where it starts)
z <- rep(1, n) # Initialize z for the GIG distribution
prior_mean_beta <- rep(0, ncol(X)+1) # Prior mean for beta
prior_prec_beta <- diag(1/50, ncol(X)+1) # Prior precision for beta

# Storage Matrices
gamma_store = rep(NA, N_iteration) # Initialize the shape parameters ~ tail heaviness
delta_store = rep(NA, N_iteration) # Initialize the scale parameter ~ spread
beta_store <- matrix(NA, nrow = N_iteration, ncol = ncol(X)+1) # Initialize the matrix to store the beta coefficients for each iteration, we have ncol(X)+1 because we have the intercept and the slope
z_store <- matrix(NA,nrow=N_iteration, ncol = n) # Store z values for each iteration


# Priors
omega <- 1 # Prior how strongly we believe in eta 20-40 balanced, 100 sticks more to prior, 5 trusts the data prior guess for variance
eta <- var(y) # average variance so we use var(y) to adapt to all stocks
xi <- 0.01 # for φ|μ Shape of the Gamma prior
x_prior <- 0.01 # for φ||μ Rate of the Gamma prior

gamma <- 1 # for μ|φ STARTING POINT OF GAMMA IN THE GIG DISTRIBUTION
delta <- 1 # for μ|φ STARTING POINT OF DELTA IN THE GIG DISTRIBUTION
phi_GAMMA <- 1 # for μ|φ STARTING POINT OF PHI IN THE GAMMA DISTRIBUTION

m_GIG <- delta/gamma # prior mean for the μ parameter in the IG distribution
nu <- n + 2*xi # for φ||μ
u2 <- n + omega - x_prior # you kfor φ|μ


#diagnosis
mu_ig_store <- rep(NA, N_iteration)
phi_store <- rep(NA, N_iteration)

for (i in 1:N_iteration){
  
  # Z draws
  residuals <- y - X %*% b[1:ncol(X)]
  q_val <- 1 + (residuals/delta)^2
  alpha_param <- sqrt(gamma^2 + b[ncol(X)+1]^2)
  
  z <- numeric(n)
  for(j in 1:n) {
    z[j] <- rgig(1, lambda = -1, 
                 chi = delta^2 * q_val[j], 
                 psi = alpha_param^2)
  }
  
  # d,D matrices and coefficient generation
  A1 <- cbind(X, z)
  w <- as.vector(1/z)
  A1w <- A1 * w
  C2inv <- prior_prec_beta 
  
  D_inv <- t(A1w) %*% A1 + C2inv
  D <- solve(D_inv)
  d <- t(A1) %*% (y*w) + C2inv %*% prior_mean_beta
  
  b <- mvrnorm(1, mu = D %*% d, Sigma = D)
  
  # Use scheme 1 from paper the Gamma-GIG
  u1 <- sum(z) + omega*eta 
  u3 <- sum(1/z) + omega/eta
  
  m_GIG <- rgig(1, lambda = (n-1)/2, chi = phi_GAMMA*u1, psi = phi_GAMMA*u3)
  phi_GAMMA <- rgamma(1,shape=(nu+1)/2,rate=u1/(2*m_GIG)-u2+u3*m_GIG/2) 
  
  #Convert to Gamma,Delta
  gamma <- sqrt(phi_GAMMA / m_GIG)
  delta <- sqrt(m_GIG * phi_GAMMA)
  
  mu_ig_store[i] <- m_GIG
  phi_store[i] <- phi_GAMMA
  
  
  beta_store[i,] <- b
  gamma_store[i] <- gamma
  delta_store[i] <- delta
  z_store[i,] <- z
  
}


# Discard burn-in samples and thin the chain
keep_idx <- seq(burnin + thin, N_iteration, by = thin)
beta_post <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]
z_post_mean <- colMeans(z_store[keep_idx, ])

# Posterior Distribution

predictions <- numeric(n_pred)
beta_median <- apply(beta_post, 2, median)# Post median
gamma_median <- median(gamma_post)# Post median
delta_median <- median(delta_post)# Post median
alpha_median <- sqrt(gamma_median^2 + beta_median[2]^2)

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

hist(y, breaks = 100, probability = TRUE,
     col = rgb(0.7, 0.7, 0.7, 0.6), border = "white",
     main = paste(ticker, "Daily Returns: Data vs Fitted NIG Posterior"),
     xlab = "Daily Return (%)",
     ylab = "Density")

# Add empirical density
lines(density(y, adjust = 1.5), col = "black", lwd = 3)

# Add fitted NIG density
x_seq <- seq(min(y) - 5, max(y) + 5, length.out = 1000)
nig_density <- dnig(x_seq,
                    mu = beta_median[1],
                    beta = beta_median[2],
                    delta = delta_median,
                    alpha = alpha_median)

lines(x_seq, nig_density, col = "red", lwd = 3, lty = 2)

legend("topright",
       legend = c("Empirical Data", "Fitted NIG"),
       col = c("black", "red"),
       lwd = c(3, 3),
       lty = c(1, 2),
       bty = "n",
       cex = 1.2)


# Predictive 

for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5,
                chi = delta_median^2,
                psi = gamma_median^2 + beta_median[2]^2)
  mu_j <- beta_median[1] + beta_median[2] * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}

cat("Real TSLA:\n")
cat("  SD:        ", sd(y), "%\n")
cat("  Skewness:  ", skewness(y), "\n")
cat("  Kurtosis:  ", kurtosis(y), "\n\n")

cat("Predicted:\n")
cat("  SD:        ", sd(predictions), "%\n")
cat("  Skewness:  ", skewness(predictions), "\n")
cat("  Kurtosis:  ", kurtosis(predictions), "\n")

# Predictive Plot

par(mfrow = c(1,1), mar = c(5,5,4,2))
hist(y, breaks = 100, probability = TRUE,
     col = rgb(0.7, 0.7, 0.7, 0.6), border = "white",
     main = paste(ticker, "Returns: Real vs Predicted"),
     xlab = "Daily Return (%)",
     xlim = c(min(c(y, predictions)), max(c(y, predictions))))
hist(predictions, breaks = 100, probability = TRUE,
     col = rgb(1, 0, 0, 0.3), border = NA, add = TRUE)
lines(density(y, adjust = 1.5), col = "black", lwd = 3)
lines(density(predictions, adjust = 1.5), col = "red", lwd = 3, lty = 2)
legend("topright", 
       legend = c("Real Data", "Predicted"),
       col = c("black", "red"), 
       lwd = c(3, 3), lty = c(1, 2),
       bty = "n")

# VaR calculation of the predictive distribution at 95% confidence level
VaR_5 <- quantile(predictions, probs = 0.05) # left side downside risk
VaR_95 <- quantile(predictions, probs = 0.95) # right side upside risk
VaR_1 <- quantile(predictions, probs = 0.01) # left side extreme downside risk
ES_5 <- mean(predictions[predictions <= VaR_5]) # Expected Shortfall at 5% 

# Outputs
cat("\n=== RISK METRICS ===\n")
cat("Expected Return:     ", round(mean(predictions), 3), "%\n")
cat("Volatility (SD):     ", round(sd(predictions), 3), "%\n\n")

cat("5% VaR (1 day):      ", round(VaR_5, 3), "%\n")
cat("1% VaR (1 day):      ", round(VaR_1, 3), "%\n")
cat("95% VaR (upside):    ", round(VaR_95, 3), "%\n\n")

cat("5% Expected Shortfall:", round(ES_5, 3), "%\n")


# Convergance


