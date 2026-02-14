# Estimating NIG distribution parameters from Bayesian viewpoint using MCMC Gibbs
library(quantmod)
library(GIGrvg)
library(fBasics)

# Get SPY data for the MCMC
getSymbols("SPY", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
SPY_adj <- Ad(SPY)
R_t <- (SPY_adj / lag(SPY_adj)) - 1
y <- R_t[-1] # Remove NA

# Priors
a0 <- 1 # sample size || MUST BE POSITIVE
a1 <- 8 # skewness
a2 <- 7 # mean
a3 <- 6 # volatility || MUST BE POSITIVE
a4 <- 15 # tail heaviness || MUST BE POSITIVE
mu <- mean(y) # best guess
delta <- sd(y) # best guess
beta <- 0 # neutral symmetric distribution
gamma <- 1
alpha <- sqrt(beta^2 + gamma^2) # page 3 

# MCMC Gibbs
mcmc_samples <- matrix(NA, nrow = 20000, ncol = 4) # create a matrix to store all the mcmc samples
colnames(mcmc_samples) <- c("m", "b", "d", "a") # give names to the matrix

for (i in 1:20000) {
  
  scale_val <- (y - mu)^2 + delta^2 # z parameter page 5
  psi_val <- alpha^2 # z parameter page 5
  z <- rgig(n = length(y), lambda = -1, chi = scale_val, psi = psi_val) # compute z for each y using GIG distribution page 5
  
  n_obs <- length(y) # How many observations we got
  a0_prime <- a0 + n_obs # page 3
  a1_prime <- a1 + sum(y) # page 3          
  a2_prime <- a2 + sum(y / z) # page 3      
  a3_prime <- a3 + sum(z) / 2  # page 3     
  a4_prime <- a4 + sum(1 / z) / 2 # page 3
  
  rho <- -a0_prime / (2 * sqrt(a3_prime * a4_prime)) #page 4
  sd_mu <- sqrt(1 / (2 * (1 - rho^2) * a4_prime)) # page 4
  sd_beta <- sqrt(1 / (2 * (1 - rho^2) * a3_prime)) # page 4
  mu_tilde <- (a2_prime - (a0_prime * a1_prime) / (2 * a3_prime)) * (sd_mu^2) # page 4
  beta_tilde <- (a1_prime - (a0_prime * a2_prime) / (2 * a4_prime)) * (sd_beta^2) # page 4
  
  # We know from the paper μ β have a joint posterior distribution that is Bivariate normal
  # So it means they aren't indepented thus we have ρ so we need to use cholesky transformation to get the correlated samples of μ and β
  # So we use z standard to generate 2 different normal random variables and then we use the cholesky transformation to get the correlated samples of μ and β
  # mu scales the first raw number to have the correct mean and standard deviation
  #then beta scales the second raw number to have the correct mean and standard deviation,
  #but also adds a term that is correlated with the first raw number to introduce the correlation between mu and beta
  
  z_standard <- rnorm(2)
  mu <- mu_tilde + sd_mu * z_standard[1]
  beta <- beta_tilde + sd_beta * (rho * z_standard[1] + sqrt(1 - rho^2) * z_standard[2]) # cholesky transformation
  
  d2_shape <- (a0_prime + 1) / 2 # page 4 gamma parameter
  d2_rate <- a4_prime - (a0_prime^2 / (4 * a3_prime)) # page 4 gamma parameter
  delta_sq <- rgamma(1, shape = d2_shape, rate = d2_rate) # page 4 generate for delta^2
  delta <- sqrt(delta_sq) # page 4  get the simple delta for γ|δ
  
  gamma_mean <- (a0_prime * delta) / (2 * a3_prime) # page 4 
  gamma_sd <- sqrt(1 / (2 * a3_prime)) # page 4
  
  # Simple rejection sampling for truncation:
  
  repeat { # we use repeat cause we need γ > 0
    gamma <- rnorm(1, mean = gamma_mean, sd = gamma_sd)
    if (gamma > 0) break
  }
  alpha <- sqrt(gamma^2 + beta^2)
  mcmc_samples[i, ] <- c(mu, beta, delta, alpha)
}

final_estimates <- mcmc_samples[seq(2001, 10000, by = 50), ] # thinning every 50th and burn in is 2000
colMeans(final_estimates)

# plot 

# Expanding xlim slightly beyond min/max helps see the tail behavior
x_range <- seq(min(y) - 0.02, max(y) + 0.02, length.out = 500) 
y_max <- max(hist(y, breaks = 50, plot = FALSE)$density) * 1.1

# 2. Enhanced Plotting
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
hist(y, breaks = 50, probability = TRUE, 
     main = "SPY Daily Simple Returns vs. Fitted NIG",
     xlab = "Daily Simple Return (decimal)", 
     ylab = "Density",
     col = "lightgray", border = "white",
     xlim = c(min(x_range), max(x_range)),
     ylim = c(0, y_max))

# 3. High-Resolution NIG Curve
# Using 500 points (x_range) makes the sharp peak look smoother
lines(x_range, dnig(x_range, 
                    alpha = p_means["a"], 
                    beta = p_means["b"], 
                    delta = p_means["d"], 
                    mu = p_means["m"]), 
      col = "red", lwd = 2.5)

# 4. Add the 'Identity' line (0% return) and Expected Return
abline(v = 0, col = "black", lty = 3) # Break-even line
p_gamma <- sqrt(p_means["a"]^2 - p_means["b"]^2)
expected_y <- p_means["m"] + (p_means["d"] * p_means["b"]) / p_gamma
abline(v = expected_y, col = "blue", lwd = 2, lty = 2) # Mean return

legend("topright", 
       legend = c("SPY Returns", "Fitted NIG", "Mean Return", "0% Return"), 
       col = c("lightgray", "red", "blue", "black"), 
       lwd = c(8, 2.5, 2, 1), 
       lty = c(1, 1, 2, 3),
       cex = 0.8)
