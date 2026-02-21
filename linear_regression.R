library(quantmod)
library(GIGrvg)   # rgig
library(MASS)     # mvrnorm

set.seed(42)

# ====================== DATA (SPY daily returns) ======================
getSymbols("TSLA", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
SPY_adj <- Ad(TSLA)
y <- as.numeric(dailyReturn(SPY_adj)[-1])   # simple returns, drop first NA
n <- length(y)

cat("Number of observations:", n, "\n")

# ====================== PRIORS (Gamma-GIG - same as paper) ======================
xi    <- 0.01      # phi ~ Gamma(xi, chi)
chi   <- 0.01
eta   <- 5         # prior mean for mu_ig
omega <- 0.001

# Diffuse prior on (μ, β)  -- variances = 50 as in the FTARET example
prior_mean_beta <- c(0, 0)
prior_prec_beta <- diag(1/50, 2)          # 2 parameters: μ and β_skew

# ====================== MCMC SETTINGS ======================
n_iter <- 100000
burn   <- 20000
thin   <- 50
keep_idx <- seq(burn + 1, n_iter, by = thin)

# Storage
beta_store  <- matrix(NA, n_iter, 2)      # columns: mu, beta_skew
gamma_store <- numeric(n_iter)
delta_store <- numeric(n_iter)

# Initial values
mu    <- mean(y)
beta  <- 0
gamma <- 1
delta <- 1
mu_ig <- delta / gamma
phi   <- 1

# ====================== MCMC LOOP (Gamma-GIG) ======================
for (i in 1:n_iter) {
  
  # ----- 1. Sample latent z_i -----
  residuals <- y - mu                       # y - μ   (no covariates)
  chi_z     <- residuals^2 + delta^2        # δ_GIG²
  psi_z     <- gamma^2 + beta^2             # α²
  z <- rgig(n = n, lambda = -1, chi = chi_z, psi = psi_z)
  
  # ----- 2. Update (μ, β) via heteroscedastic regression -----
  X_full <- cbind(1, z)                     # design matrix: intercept + z
  w      <- 1 / z
  D_inv  <- t(X_full * w) %*% X_full + prior_prec_beta
  D      <- solve(D_inv)
  d      <- t(X_full) %*% (y * w) + prior_prec_beta %*% prior_mean_beta
  b      <- mvrnorm(1, mu = D %*% d, Sigma = D)
  
  mu   <- b[1]
  beta <- b[2]
  
  # ----- 3. Gamma-GIG update for IG parameters -----
  z_bar   <- mean(z)
  z_bar_r <- mean(1/z)
  nu      <- n + 2 * xi
  u1      <- n * z_bar + omega * eta
  u2      <- n + omega - chi
  u3      <- n * z_bar_r + omega / eta
  
  mu_ig <- rgig(1, lambda = (n-1)/2, chi = phi * u1, psi = phi * u3)
  
  rate_phi <- u1/(2*mu_ig) - u2 + u3*mu_ig/2
  phi      <- rgamma(1, shape = (nu + 1)/2, rate = rate_phi)
  
  # Back-transform to usual NIG parameters
  gamma <- sqrt(phi / mu_ig)
  delta <- sqrt(mu_ig * phi)
  
  # Store
  beta_store[i, ] <- c(mu, beta)
  gamma_store[i]  <- gamma
  delta_store[i]  <- delta
}

# ====================== POSTERIOR SUMMARY ======================
beta_post  <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]

cat("\n=== Posterior medians (Gamma-GIG) ===\n")
cat("μ     =", round(median(beta_post[,1]), 5), "\n")
cat("β     =", round(median(beta_post[,2]), 5), "\n")
cat("γ     =", round(median(gamma_post), 5), "\n")
cat("δ     =", round(median(delta_post), 5), "\n")
cat("α     =", round(median(sqrt(gamma_post^2 + beta_post[,2]^2)), 5), "\n")

# ====================== ONE-DAY-AHEAD PREDICTIVE DISTRIBUTION ======================
n_pred <- 10000
predictions <- numeric(n_pred)

for (j in 1:n_pred) {
  idx  <- sample(1:length(keep_idx), 1)
  mu_j <- beta_post[idx, 1]
  beta_j <- beta_post[idx, 2]
  gamma_j <- gamma_post[idx]
  delta_j <- delta_post[idx]
  
  # New mixing variable from IG(γ, δ)
  z_new <- rgig(1, lambda = -0.5, chi = delta_j^2, psi = gamma_j^2)
  
  # Predictive draw
  mu_new <- mu_j + beta_j * z_new
  predictions[j] <- rnorm(1, mean = mu_new, sd = sqrt(z_new))
}

cat("\nTomorrow's SPY return forecast (posterior predictive):\n")
cat("Expected return :", round(mean(predictions)*100, 3), "%\n")
cat("95% credible interval :", round(quantile(predictions, c(0.025, 0.975))*100, 3), "%\n")
cat("5% VaR (loss) :", round(quantile(predictions, 0.05)*100, 3), "%\n")