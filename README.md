# Bayesian NIG Model for Stock Returns

A simple implementation of the Normal-Inverse Gaussian (NIG) distribution using Bayesian MCMC methods to model stock returns and estimate Value-at-Risk (VaR).

---

## 📌 What Does This Do?

This code fits a **Normal-Inverse Gaussian (NIG) distribution** to stock returns using **Bayesian MCMC (Gibbs sampling)**. 

**Why?** Because normal distributions are **terrible** at predicting tail risk in financial markets. NIG captures fat tails and skewness much better.

---

## 🎯 Main Results

**Tested on S&P 500 Total Return (2015-2025):**

| Metric | Real Data | NIG Model | Normal Distribution |
|--------|-----------|-----------|---------------------|
| **1% VaR** | -3.32% | -3.38% ✅ | -2.60% ❌ |
| **5% VaR** | -1.88% | -1.78% ✅ | -1.58% ❌ |
| **Kurtosis** | 15.51 | 8.51 ✅ | 3.00 ❌ |

**Bottom line:** NIG gives much more realistic tail risk estimates than Normal distribution.

---

## 📐 The Model

### Simple Version

```
Daily Return = μ + β·z + √z·ε

where:
  μ = average return
  β = skewness parameter
  z = latent variance (different each day)
  ε = random noise
```

### Why This Works

- **Normal distribution:** Assumes same variance every day (unrealistic)
- **NIG model:** Each day gets its own variance `z` (realistic!)
  - Calm days: small z
  - Crisis days: large z

---

## 🔧 How It Works

### 1. Data Preparation
```r
# Get S&P 500 data
ticker <- "^SP500TR"
data <- getSymbols(ticker, from = "2015-12-01")
y <- 100 * dailyReturn(data)  # Convert to percent
```

### 2. MCMC Estimation
```r
# Run Gibbs sampler
For 500,000 iterations:
  1. Sample z (variance for each day)
  2. Sample μ, β (regression coefficients)
  3. Sample γ, δ (NIG parameters)
```

### 3. Generate Predictions
```r
# Simulate future returns
For 50,000 predictions:
  - Draw random z from fitted distribution
  - Generate return using fitted parameters
```

### 4. Calculate VaR
```r
VaR_1 <- quantile(predictions, 0.01)  # 1% worst case
VaR_5 <- quantile(predictions, 0.05)  # 5% worst case
```

---

## 📊 What You Get

### Parameter Estimates
```
γ (gamma) = 0.577   — Controls tail heaviness (smaller = fatter tails)
δ (delta) = 0.674%  — Volatility scale
μ (mu)    = 0.056%  — Average return
β (beta)  = -0.020  — Skewness (negative = left-skewed)
```

### Risk Metrics
```
Expected Return:      0.06%
Volatility (SD):      1.06%

5% VaR (1 day):      -1.78%  → 5% chance of losing more than this
1% VaR (1 day):      -3.38%  → 1% chance of losing more than this
Expected Shortfall:  -2.75%  → Average loss in worst 5% scenarios
```

---

## 🚀 How to Use

### Requirements
```r
install.packages(c("quantmod", "GIGrvg", "MASS", "fBasics"))
```

### Run the Code
```r
# 1. Load the script
source("NIG_Code_hetero.R")

# 2. The code automatically:
#    - Downloads S&P 500 data
#    - Runs MCMC estimation
#    - Generates predictions
#    - Calculates VaR
#    - Creates plots

# 3. Check results in console output
```

### Change the Stock
```r
# Edit this line in garc.R:
ticker <- "SPY"    # Or "AAPL", "TSLA", etc.
```

---

## 📈 Files in This Repo

- **NIG_Code_hetero.R** — Main code (MCMC estimation + VaR calculation)
- **NIG_Normal.R** — Normal NIG
- **README.md** — This file

---

## 🔍 Key Settings

```r
N_iteration <- 500000   # MCMC iterations
burnin <- 150000        # Discard first 150k (warm-up)
thin <- 200             # Keep every 200th sample
n_pred <- 50000         # Number of predictions

omega <- 1              # Prior strength (weak = trust data)
eta <- var(y)           # Prior scale (adaptive)
```

**Runtime:** ~20 minutes on modern laptop

---

## 📚 What Is This Based On?

**Paper:** Karlis & Lillestøl (2004) - "Bayesian estimation of NIG models via Markov chain Monte Carlo methods"

**Method:** Gibbs sampling with Gamma-GIG scheme

**Innovation:** Heteroscedastic regression (each day has different variance)

---

## ✅ Advantages Over Normal Distribution

| Feature | Normal | NIG |
|---------|--------|-----|
| Tail risk | ❌ Underestimates | ✅ Realistic |
| Kurtosis | ❌ Always 3 | ✅ Flexible (matches data) |
| Skewness | ❌ Always 0 | ✅ Captures asymmetry |
| Crisis periods | ❌ Misses | ✅ Adapts via latent variance |

---

## 🎯 Practical Use Cases

1. **Risk Management:** Better VaR estimates for portfolios
2. **Regulatory Capital:** More realistic worst-case scenarios
3. **Option Pricing:** Fat-tailed returns affect option values
4. **Portfolio Optimization:** Account for realistic tail risk

---

## 🔮 Future Extensions

- Add covariates (VIX, interest rates) for conditional forecasting
- Test on other asset classes (crypto, commodities, FX)
- Out-of-sample backtesting
- Compare with GARCH models

---

## 📧 Contact

**Author:** Scaply  
**University:** Piraeus - Finance & Banking Management  
**GitHub:** [github.com/Scuplex](https://github.com/Scuplex)

Questions? Open an issue!

---

## 📝 License

MIT License - Free for academic and research use

---

**Last Updated:** March 2026
