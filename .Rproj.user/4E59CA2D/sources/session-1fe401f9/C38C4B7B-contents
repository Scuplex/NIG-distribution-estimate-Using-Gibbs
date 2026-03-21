# =============================================
# FINAL POLISHED VERSION - NIG + AR(1) + Proper Backtest
# =============================================
library(quantmod)
library(GIGrvg)
library(MASS)
library(fBasics)

# ==================== DATA ====================
ticker <- "SPY"   # Change to "TSLA" anytime

data_raw <- getSymbols(ticker, from = "2015-12-01", to = "2025-12-31", 
                       periodicity = "daily", auto.assign = FALSE)

ret <- 100 * as.numeric(dailyReturn(Ad(data_raw), type = "arithmetic"))
y_full <- ret[-1]
lag_full <- ret[-length(ret)]

# Train/Test split (80/20)
train_size <- floor(0.8 * length(y_full))
y_train   <- y_full[1:train_size]
lag_train <- lag_full[1:train_size]
y_test    <- y_full[(train_size+1):length(y_full)]
lag_test  <- lag_full[(train_size+1):length(y_full)]

n <- length(y_train)

# ==================== SETTINGS ====================
N_iteration <- 80000
burnin      <- 20000
thin        <- 80

# Design matrix: intercept + AR(1)
X <- cbind(1, lag_train)
b <- rep(0, ncol(X) + 1)

# Universal Priors
xi        <- 0.01
x_prior   <- 0.01
phi_GAMMA <- 1
eta       <- var(y_train) * 1.0
omega     <- 30          # slightly stronger for stability
gamma     <- 1
delta     <- 1
m_GIG     <- 1

nu <- n + 2*xi
u2 <- n + omega - x_prior

# Storage
beta_store  <- matrix(NA, N_iteration, ncol(X)+1)
gamma_store <- delta_store <- rep(NA, N_iteration)

cat("Training AR(1) NIG model on", ticker, "...\n")

# ==================== MCMC LOOP ====================
for (i in 1:N_iteration) {
  
  # Z draws - Correct & Fast
  residuals   <- y_train - X %*% b[1:ncol(X)]
  q_val       <- 1 + (residuals / delta)^2
  alpha_param <- sqrt(gamma^2 + b[ncol(X)+1]^2)
  
  z <- sapply(q_val, function(q) rgig(1, -1, delta^2 * q, alpha_param^2))
  
  # Regression update (AR(1) + skew)
  A1 <- cbind(X, z)
  w  <- 1/z
  C2inv <- diag(1/50, ncol(A1))
  D_inv <- t(A1 * w) %*% A1 + C2inv
  D <- solve(D_inv)
  d <- t(A1) %*% (y_train * w) + C2inv %*% rep(0, ncol(A1))
  b <- mvrnorm(1, D %*% d, D)
  
  # Gamma-GIG update
  u1 <- sum(z) + omega * eta
  u3 <- sum(1/z) + omega / eta
  m_GIG <- rgig(1, (n-1)/2, phi_GAMMA*u1, phi_GAMMA*u3)
  phi_GAMMA <- rgamma(1, (nu+1)/2, rate = u1/(2*m_GIG) - u2 + u3*m_GIG/2)
  
  gamma <- sqrt(phi_GAMMA / m_GIG)
  delta <- sqrt(m_GIG * phi_GAMMA)
  
  beta_store[i,]  <- b
  gamma_store[i]  <- gamma
  delta_store[i]  <- delta
}

# Extract posterior
keep_idx <- seq(burnin + thin, N_iteration, by = thin)
beta_post  <- beta_store[keep_idx, ]
gamma_post <- gamma_store[keep_idx]
delta_post <- delta_store[keep_idx]

beta_median <- apply(beta_post, 2, median)
gamma_median <- median(gamma_post)
delta_median <- median(delta_post)

cat("\n=== TRAINED PARAMETERS ===\n")
cat("Intercept     :", round(beta_median[1], 5), "\n")
cat("AR(1) coef    :", round(beta_median[2], 5), "\n")
cat("Skewness beta :", round(beta_median[3], 5), "\n")
cat("gamma         :", round(gamma_median, 5), "\n")
cat("delta         :", round(delta_median, 5), "\n\n")

# ==================== PROPER OUT-OF-SAMPLE PREDICTION ====================
cat("Generating proper out-of-sample predictions...\n")

n_pred_per_day <- 5000
predictions_matrix <- matrix(NA, length(y_test), 5)
colnames(predictions_matrix) <- c("Real", "Pred_Mean", "Pred_SD", "VaR5", "VaR95")

current_mu_forecast <- beta_median[1] + beta_median[2] * lag_test[1]

for (t in 1:length(y_test)) {
  preds <- numeric(n_pred_per_day)
  
  for (j in 1:n_pred_per_day) {
    q_i <- 1 + ((current_mu_forecast - beta_median[1] - beta_median[2]*lag_test[t]) / delta_median)^2
    
    z_new <- rgig(1, -1, delta_median^2 * q_i, gamma_median^2 + beta_median[3]^2)
    
    mu_j <- beta_median[1] + beta_median[2]*lag_test[t] + beta_median[3] * z_new
    preds[j] <- rnorm(1, mu_j, sqrt(z_new))
  }
  
  predictions_matrix[t, "Real"]      <- y_test[t]
  predictions_matrix[t, "Pred_Mean"] <- mean(preds)
  predictions_matrix[t, "Pred_SD"]   <- sd(preds)
  predictions_matrix[t, "VaR5"]      <- quantile(preds, 0.05)
  predictions_matrix[t, "VaR95"]     <- quantile(preds, 0.95)
  
  current_mu_forecast <- mean(preds)   # update for next day
}

# ==================== RESULTS ====================
cat("\n=== OUT-OF-SAMPLE PERFORMANCE ===\n")
cat("MSE  :", round(mean((predictions_matrix[,"Real"] - predictions_matrix[,"Pred_Mean"])^2), 4), "\n")
cat("MAE  :", round(mean(abs(predictions_matrix[,"Real"] - predictions_matrix[,"Pred_Mean"])), 4), "\n")

viol5  <- mean(predictions_matrix[,"Real"] < predictions_matrix[,"VaR5"])
viol95 <- mean(predictions_matrix[,"Real"] > predictions_matrix[,"VaR95"])
cat("5% VaR violation rate  :", round(100*viol5, 1), "%\n")
cat("95% VaR violation rate :", round(100*viol95, 1), "%\n")

# ==================== PLOTS (exactly like your image) ====================
par(mfrow = c(2,2), mar = c(4,4,3,2))

# Plot 1: Realized vs Predicted
plot(predictions_matrix[,"Real"], predictions_matrix[,"Pred_Mean"], 
     pch=16, col=rgb(0,0,1,0.4), main="Realized vs Predicted", 
     xlab="Realized Return (%)", ylab="Predicted Mean (%)")
abline(0,1, col="red", lwd=2)

# Plot 2: Prediction Errors
errors <- predictions_matrix[,"Real"] - predictions_matrix[,"Pred_Mean"]
plot(errors, type="l", main="Prediction Errors Over Time", ylab="Error (%)")
abline(h=0, col="red", lty=2)

# Plot 3: VaR Coverage
plot(1:length(y_test), predictions_matrix[,"Real"], pch=16, col="black", cex=0.6,
     main="VaR Coverage (Out-of-Sample)", xlab="Test Day", ylab="Return (%)")
lines(1:length(y_test), predictions_matrix[,"VaR5"], col="red", lwd=2)
lines(1:length(y_test), predictions_matrix[,"VaR95"], col="green", lwd=2)
legend("topright", c("Realized", "5% VaR", "95% VaR"), 
       col=c("black","red","green"), lwd=c(1,2,2), pch=c(16,NA,NA))

# Plot 4: Density Comparison
hist(y_test, breaks=60, probability=TRUE, col=rgb(0.7,0.7,0.7,0.6), 
     main="Test Set: Real vs Model Distribution", xlab="Return (%)")
lines(density(predictions_matrix[,"Pred_Mean"], adjust=1.5), col="red", lwd=3, lty=2)
legend("topright", c("Real Test Data", "Model Predictions"), 
       col=c("gray","red"), lwd=c(10,3), lty=c(1,2))

cat("\n✅ DONE! Your final AR(1) NIG backtest is ready.\n")