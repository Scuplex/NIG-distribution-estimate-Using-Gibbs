# Bayesian estimation of heteroscedastic regression

library(quantmod)
library(GIGrvg)
library(MASS)
library(parallel)
library(fBasics)
#set.seed(42) # for reproducibility

# Get SPY data for the MCMC

ticker <- "QQQ"
start_date <- "2015-12-01"
end_date <- "2025-12-31"

data_raw <- getSymbols(ticker, from = start_date, to = end_date, periodicity = "daily", auto.assign = FALSE)
y <- 100 * as.numeric(dailyReturn(Ad(data_raw), type = "arithmetic"))[-1]  # percent returns
n <- length(y)

N_iteration <- 300000
burnin <- 40000
thin  <- 120
n_pred <- 50000

# Covariates


# Initialize Variables/Matrices
X <- matrix(1, nrow = n, ncol = 1) # Add a column of ones for the intercept for now we only have 1 for intercept
b <- rep(0, ncol(X)+1) # Initialize the coefficients vector (0 cause we don't know where it starts)
z <- rep(1, n) # Initialize z for the GIG distribution
prior_mean_beta <- rep(0, ncol(X)+1) # Prior mean for beta
prior_prec_beta <- diag(1/50, ncol(X)+1) # Prior precision for beta

# Storage Matrices
gamma_store = rep(NA, N_iteration) # Initialize the shape parameters
delta_store = rep(NA, N_iteration) # Initialize the scale parameter
beta_store <- matrix(NA, nrow = N_iteration, ncol = ncol(X)+1) # Initialize the matrix to store the beta coefficients for each iteration, we have ncol(X)+1 because we have the intercept and the slope
z_store <- matrix(NA,nrow=N_iteration, ncol = n) # Store z values for each iteration


# Priors
omega <- 30 # Prior how strongly we believe in eta 20-40 balanced, 100 sticks more to prior, 5 trusts the data
eta <- var(y) * 1.0 # average variance so we use var(y) to adapt to all stocks
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
  
  beta_store[i,] <- b
  gamma_store[i] <- gamma
  delta_store[i] <- delta
  z_store[i,] <- z
  
}


# Discard and Thinning
keep_idx <- seq(burnin + thin, N_iteration, by = thin)
beta_post <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]
z_post_mean <- colMeans(z_store[keep_idx, ])


# Posterior Distribution

predictions <- numeric(n_pred) # Post median
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



# Predictive Distribution

for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5,
                chi = delta_median^2,
                psi = gamma_median^2 + beta_median[2]^2)
  mu_j <- beta_median[1] + beta_median[2] * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}

# POSTERIOR PREDICTIVE + VaR
cat("\nPOSTERIOR PREDICTIVE + VaR\n")
cat("Statistic| Real Data | Model Predictions\n")
cat(sprintf("Mean                | %8.4f    | %8.4f\n", mean(y), mean(predictions)))
cat(sprintf("SD                  | %8.4f    | %8.4f\n", sd(y), sd(predictions)))
cat(sprintf("Skewness            | %8.4f    | %8.4f\n", skewness(y), skewness(predictions)))
cat(sprintf("Excess Kurtosis     | %8.4f    | %8.4f\n", kurtosis(y)-3, kurtosis(predictions)-3))
cat(sprintf("5%% VaR              | %8.4f    | %8.4f\n", quantile(y, 0.05), quantile(predictions, 0.05)))
cat(sprintf("95%% VaR             | %8.4f    | %8.4f\n", quantile(y, 0.95), quantile(predictions, 0.95)))

# viol_rate <- mean(y < quantile(predictions, 0.05))
# cat(sprintf("5%% VaR Violation Rate | %.2f%%     (ideal ≈ 5%%)\n", viol_rate*100))

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
