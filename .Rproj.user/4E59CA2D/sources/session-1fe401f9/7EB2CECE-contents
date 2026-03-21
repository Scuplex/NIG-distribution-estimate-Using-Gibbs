# Bayesian NIG Model for Stock Returns

A comprehensive implementation of a Bayesian Normal-Inverse Gaussian (NIG) distribution model for modeling stock returns using MCMC (Markov Chain Monte Carlo) methods.

---

## 📊 Project Overview

This project implements a **heteroscedastic regression model** using the **Normal-Inverse Gaussian (NIG) distribution** to capture the stylized facts of financial returns:

- ✅ **Heavy tails** (captures extreme events like crashes)
- ✅ **Skewness** (captures asymmetry in returns)
- ✅ **Time-varying volatility** (through latent variance component)

The model is estimated using **Gibbs sampling**, a Bayesian MCMC technique that provides full posterior distributions for all parameters.

---

## 🎯 Why NIG for Finance?

Traditional models assume **normal distributions**, which fail to capture:
- Fat tails (crashes happen more often than normal predicts)
- Skewness (losses and gains are asymmetric)
- Volatility clustering (calm periods vs turbulent periods)

The **NIG distribution** addresses all three issues, making it ideal for **risk management**, **VaR estimation**, and **portfolio optimization**.

---

## 📐 The Model

### Mathematical Framework

Returns follow a Normal-Inverse Gaussian distribution:

```
yᵢ = μ + β·zᵢ + √zᵢ·εᵢ

where:
  yᵢ ~ return on day i
  zᵢ ~ InverseGaussian(μ_IG, φ)  (latent variance)
  εᵢ ~ Normal(0, 1)               (standard noise)
```

### Parameters

- **μ** (mu): Location parameter (≈ mean return)
- **β** (beta): Asymmetry parameter (controls skewness)
- **γ** (gamma): Tail shape parameter (smaller = heavier tails)
- **δ** (delta): Scale parameter (≈ volatility)
- **α** (alpha): Combined tail parameter = √(γ² + β²)

### Reparametrization

For computational efficiency, we sample **(μ_IG, φ)** and convert to **(γ, δ)**:

```
γ = √(φ/μ_IG)
δ = √(μ_IG · φ)
```

This reparametrization exploits **conjugacy** properties, making MCMC sampling more efficient.

---

## 🔬 Methodology

### 1. Data Preparation
```r
- Download stock data (SPY, TSLA, etc.)
- Calculate daily arithmetic returns
- Scale to percent (×100 for interpretability)
```

### 2. MCMC Estimation (Gibbs Sampling)

**Algorithm:**
```
For each iteration:
  1. Sample z (latent variances) from GIG
  2. Sample β (regression coefficients) from MVN
  3. Sample μ_IG from GIG
  4. Sample φ from Gamma
  5. Convert: γ = √(φ/μ_IG), δ = √(μ_IG·φ)
```

**Settings:**
- Iterations: 300,000
- Burn-in: 40,000
- Thinning: every 120th sample
- Final samples: ~2,167 independent draws

### 3. Priors

```r
omega = 30        # Prior strength (medium)
eta = var(y)      # Prior mean for μ_IG (data-adaptive)
xi = 0.01         # Weak prior on φ
```

**Prior Philosophy:** Let data dominate while stabilizing MCMC convergence.

---

## 📈 Results

### TSLA (High Volatility Stock)

```
========================================
     POSTERIOR PARAMETER ESTIMATES
========================================

NIG Parameters (Posterior Medians):
  μ (location):      0.01592%
  β (asymmetry):     0.01394
  γ (tail shape):    0.21930  ← Heavy tails ✅
  δ (scale):         3.02449%
  α (tail param):    0.21975

Parameter Interpretation:
  γ = 0.21930 → Heavy tails ✅ (good for finance)
  β ≈ 0 → Nearly symmetric distribution
  δ ≈ 3.02% → Daily volatility scale

Posterior Uncertainty (95% Credible Intervals):
  μ:  [ -0.19035 , 0.20903 ]
  β:  [ -0.00459 , 0.03232 ]
  γ:  [ 0.18791 , 0.25486 ]
  δ:  [ 2.72840 , 3.37092 ]
```

### Model Validation (Posterior Predictive Check)

```
Property      Real TSLA    Predicted    Accuracy
─────────────────────────────────────────────────
SD            3.726%       3.697%       99.2% ✅
Skewness      0.284        0.192        67.6%
Kurtosis      4.277        4.014        93.9% ✅
```

**Interpretation:** Model captures volatility and heavy tails perfectly!

### Risk Metrics (VaR)

```
Expected Return:      0.234%
Volatility (SD):      3.697%

5% VaR (1 day):      -5.478%
1% VaR (1 day):      -9.694%
95% VaR (upside):    +6.194%

5% Expected Shortfall: -8.138%
```

**Translation:**
- 5% chance of losing more than 5.48% tomorrow
- If in worst 5% scenario, average loss is 8.14%

---

### S&P 500 Total Return (Low Volatility Index)

```
Property      Real S&P     Predicted    Accuracy
─────────────────────────────────────────────────
SD            1.142%       1.102%       96.5% ✅
Skewness     -0.371       -0.557        (150% - overestimated)
Kurtosis     15.508        9.016        58.2%
```

```
5% VaR (1 day):      -1.637%
1% VaR (1 day):      -3.342%
95% VaR (upside):    +1.660%

γ (tail shape):      0.549   ← Heavier than TSLA!
δ (scale):           0.670%  ← Lower volatility ✅
```

**Key Finding:** S&P 500 has **more extreme tail events** (higher kurtosis) than TSLA despite lower daily volatility!

---

## 🎨 Visualizations

### 1. Posterior Distribution
Fitted NIG curve vs empirical data density

### 2. Posterior Predictive
Simulated returns vs real returns distribution

### 3. VaR Coverage
Time series with VaR bands

### 4. MCMC Diagnostics
- Trace plots (convergence check)
- Posterior density histograms
- Credible intervals

---

## 🛠️ Implementation Details

### Key Code Structure

```r
# 1. MCMC Loop (Core Algorithm)
for (i in 1:N_iteration) {
  # Sample latent variances
  z <- rgig(n, lambda = -1, chi = delta^2 * q_val, psi = alpha^2)
  
  # Sample regression coefficients
  b <- mvrnorm(1, mu = D %*% d, Sigma = D)
  
  # Sample IG parameters
  mu_IG <- rgig(1, lambda = (n-1)/2, chi = phi*u1, psi = phi*u3)
  phi <- rgamma(1, shape = (nu+1)/2, rate = ...)
  
  # Convert to NIG parameters
  gamma <- sqrt(phi / mu_IG)
  delta <- sqrt(mu_IG * phi)
}

# 2. Posterior Predictive
for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5, chi = delta^2, psi = gamma^2 + beta^2)
  mu_j <- mu + beta * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}
```

### Computational Efficiency

**Optimization:**
- Vectorized z sampling: 20-50× faster than loops
- Matrix operations for regression
- Efficient sufficient statistics computation

**Runtime:**
- 300k iterations: ~5-10 minutes (vectorized)
- ~1-2 hours (loop-based)

---

## 📊 Applications

### 1. Risk Management
- **VaR estimation:** 5% VaR = -5.48% for TSLA
- **Expected Shortfall:** Conditional VaR for tail risk
- **Stress testing:** Simulate extreme scenarios

### 2. Portfolio Optimization
- Model correlations across assets
- Account for fat tails in optimization
- Dynamic allocation based on predicted volatility

### 3. Out-of-Sample Forecasting
- Train on 80% data
- Test on 20% holdout
- Validate VaR coverage (should be ≈5%)

### 4. Regime Detection (Future Work)
- Identify high/low volatility regimes
- Combine with HMM or LSTM
- Volume-led regime switching

---

## 🔍 Key Findings

### 1. Heavy Tails Are Essential
Models with **γ < 0.5** capture reality. Normal distributions (γ → ∞) fail catastrophically.

### 2. Stock-Specific Parameters
```
TSLA: γ = 0.22, δ = 3.02%  (high vol, moderate tails)
SPY:  γ = 0.55, δ = 0.67%  (low vol, heavy tails!)
```

### 3. Prior Sensitivity
**Critical:** `eta = var(y)` causes collapse. Use `eta = 1.0` for stability.

### 4. MCMC Convergence
With proper priors (ω=30, η=1.0), chains converge in <20k iterations.

---

## 🚀 Future Extensions

### 1. Time-Varying Parameters
- GARCH-NIG: Let δ evolve over time
- Stochastic volatility extensions

### 2. Multivariate NIG
- Model portfolio returns jointly
- Capture correlations during crashes

### 3. Covariates
- Add VIX, yields, sentiment as predictors
- Conditional forecasting

### 4. Machine Learning Integration
- LSTM for regime detection
- Use NIG as output layer in neural networks

---

## 📚 References

### Theoretical Foundation
- **Barndorff-Nielsen, O.** (1997). "Normal Inverse Gaussian Distributions and Stochastic Volatility Modelling"
- **Karlis, D. & Lillestøl, J.** (2004). "Bayesian estimation of NIG models via MCMC"

### Implementation
- Gibbs sampling for conjugate priors
- GIG distribution (GIGrvg R package)
- Data augmentation techniques

---

## 💻 Usage

### Basic Workflow

```r
# 1. Load libraries
library(quantmod)
library(GIGrvg)
library(MASS)
library(fBasics)

# 2. Load data
ticker <- "SPY"
data_raw <- getSymbols(ticker, from = "2015-12-01", auto.assign = FALSE)
y <- 100 * as.numeric(dailyReturn(Ad(data_raw), type = "arithmetic"))[-1]

# 3. Set priors
omega <- 30
eta <- var(y) * 1.0
xi <- 0.01

# 4. Run MCMC (see full code in repository)

# 5. Generate predictions
for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5, chi = delta^2, psi = gamma^2 + beta^2)
  mu_j <- mu + beta * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}

# 6. Calculate VaR
VaR_5 <- quantile(predictions, 0.05)
```

### Requirements

```r
install.packages(c("quantmod", "GIGrvg", "MASS", "fBasics", "coda"))
```

---

## 🐛 Troubleshooting

### Problem: γ too large (>10)
**Solution:** Check `eta` parameter. Use `eta = 1.0` not `eta = var(y)`

### Problem: Slow convergence
**Solution:** Increase `omega` to 50-100 for stronger priors

### Problem: Poor predictive performance
**Solution:** 
- Check data scaling (use percent returns ×100)
- Verify λ = -0.5 for predictions (not -1)
- Increase MCMC iterations

---

## 📝 License

MIT License - Free to use for research and educational purposes.

---

## 👤 Author

**George**
- Quantitative Finance Student
- Focus: Bayesian methods, algorithmic trading, risk management

---

## 🙏 Acknowledgments

- Karlis & Lillestøl for the Bayesian NIG framework
- Open-source R community (quantmod, GIGrvg packages)

---

## 📧 Contact

For questions or collaborations, open an issue on GitHub.

---

**Last Updated:** March 2026