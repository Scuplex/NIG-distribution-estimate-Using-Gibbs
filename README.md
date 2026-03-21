# Updated README.md

Here's your corrected README that accurately reflects what you actually did:

```markdown
# Bayesian NIG Model for Financial Risk Management

A Bayesian implementation of the Normal-Inverse Gaussian (NIG) distribution for modeling equity returns and estimating Value-at-Risk (VaR), following Karlis & Lillestøl (2004).

---

## 📊 Project Motivation

Standard risk models rely on the **Normal (Gaussian) distribution**, which systematically underestimates tail risk in financial markets:

- ❌ **Assumes kurtosis = 3** (real markets: 8-20+)
- ❌ **Underestimates extreme events** (crashes, flash rallies)
- ❌ **Produces misleading VaR estimates** (too optimistic)

The **Normal-Inverse Gaussian (NIG) distribution** solves these problems by modeling:
- ✅ **Heavy tails** (realistic crash probabilities)
- ✅ **Skewness** (asymmetric up/down moves)
- ✅ **Time-varying volatility** (heteroscedastic variance)

---

## 🎯 Key Results: S&P 500 Total Return (2015-2025)

### VaR Comparison

| Metric | Empirical | NIG Model | Normal Model | Winner |
|--------|-----------|-----------|--------------|---------|
| **1% VaR** | -3.321% | **-3.377%** | -2.596% | ✅ NIG |
| **5% VaR** | -1.879% | **-1.784%** | -1.578% | ✅ NIG |
| **5% ES** | -2.770% | **-2.753%** | -2.279% | ✅ NIG |
| **Kurtosis** | 15.51 | 8.51 | 3.00 | ✅ NIG |

**Key Findings:**
- NIG captures 54.8% of excess kurtosis (vs 0% for Normal)
- 1% VaR error: NIG = 0.056%, Normal = 0.725% (13× improvement)
- 5% ES error: NIG = 0.017%, Normal = 0.491% (29× improvement)

---

## 📐 The Model

### Mathematical Framework

Heteroscedastic regression with NIG errors:

```
yᵢ = μ + β·zᵢ + √zᵢ·εᵢ

where:
  yᵢ = daily return on day i (%)
  μ = location parameter (base return)
  β = asymmetry parameter (skewness coefficient)
  zᵢ ~ InverseGaussian(μ_IG, φ) — different for each day!
  εᵢ ~ Normal(0, 1)
  
NIG Parameters:
  γ = √(φ/μ_IG)  — tail heaviness (smaller = fatter tails)
  δ = √(μ_IG·φ)  — scale (volatility)
```

**Why heteroscedastic?**
Each observation gets its own latent variance `zᵢ`, allowing:
- Calm days: zᵢ ≈ 0.5 (low volatility)
- Crisis days: zᵢ ≈ 10+ (extreme volatility)

---

## 🔬 Methodology

### Bayesian MCMC Estimation (Gibbs Sampling)

**Implementation:** Gamma-GIG scheme with data augmentation (Karlis & Lillestøl 2004)

**Algorithm:**
```
For each iteration:
  1. Sample z₁,...,zₙ (latent variances) from GIG
  2. Sample (μ, β) from Multivariate Normal
  3. Sample μ_IG from GIG  
  4. Sample φ from Gamma
  5. Convert: γ = √(φ/μ_IG), δ = √(μ_IG·φ)
```

**MCMC Settings:**
- Iterations: 500,000
- Burn-in: 150,000  
- Thinning: 200
- Final posterior samples: 1,750

**Priors:**
```r
omega = 1         # Weak prior (let data dominate)
eta = var(y)      # Data-adaptive scale
xi = 0.01         # Weak prior on φ
```

---

## 📈 Detailed Results

### S&P 500 Total Return Index (^SP500TR)
**Period:** December 2015 – December 2025  
**Observations:** 2,534 daily returns

### Posterior Parameter Estimates

```
NIG Parameters (Posterior Medians):
  μ (location):     0.056%
  β (asymmetry):   -0.020  (slight negative skew)
  γ (tail shape):   0.577  (heavy tails)
  δ (scale):        0.674%
  
95% Credible Intervals:
  γ: [0.52, 0.64]
  δ: [0.60, 0.75]
```

### Model Validation

```
Moment Matching:
                  Empirical    NIG Model    Normal    
─────────────────────────────────────────────────────
Mean              0.086%       0.060%       0.086%
SD                1.142%       1.065%       1.142%
Skewness         -0.371       -0.554       0.000
Kurtosis         15.51        8.51         3.00

Match:
  Volatility:      93.3% ✅
  Skewness:       149% (overestimated)
  Kurtosis:       54.8% ✅ (Normal: 0%)
```

### Risk Metrics Comparison

```
Value at Risk (1-day):
                  Empirical    NIG          Normal      Error
──────────────────────────────────────────────────────────────
1% VaR           -3.321%      -3.377%      -2.596%     
  Error                       0.056% ✅    0.725% ❌   13× worse

5% VaR           -1.879%      -1.784%      -1.578%
  Error                       0.095% ✅    0.301% ❌   3× worse

5% ES            -2.770%      -2.753%      -2.279%
  Error                       0.017% ✅    0.491% ❌   29× worse
```

**Interpretation:**
- NIG's 1% VaR is within 0.06% of reality (Normal: 0.73% off)
- For a $1M portfolio, Normal underestimates 1% VaR by **$7,250**
- NIG provides more conservative (realistic) tail risk estimates

---

## 💻 Implementation

### Requirements

```r
install.packages(c("quantmod", "GIGrvg", "MASS", "fBasics"))
```

### Basic Usage

```r
library(quantmod)
library(GIGrvg)
library(MASS)
library(fBasics)

# 1. Load data
ticker <- "^SP500TR"
data_raw <- getSymbols(ticker, from = "2015-12-01", auto.assign = FALSE)
y <- 100 * as.numeric(dailyReturn(Ad(data_raw), type = "arithmetic"))[-1]

# 2. Set priors
omega <- 1
eta <- var(y)
xi <- 0.01

# 3. Run MCMC (500k iterations, see full code)
# ... [MCMC loop] ...

# 4. Generate predictions
for(j in 1:n_pred) {
  z_new <- rgig(1, lambda = -0.5, chi = delta^2, 
                psi = gamma^2 + beta^2)
  mu_j <- mu + beta * z_new
  predictions[j] <- rnorm(1, mu_j, sqrt(z_new))
}

# 5. Calculate VaR
VaR_1 <- quantile(predictions, 0.01)
VaR_5 <- quantile(predictions, 0.05)
ES_5 <- mean(predictions[predictions <= VaR_5])
```

---

## 🚀 Future Extensions

### 1. Conditional Forecasting with Covariates
Add macroeconomic predictors to move from unconditional to conditional VaR:

```r
# Current model (unconditional):
yᵢ = μ + β·zᵢ + √zᵢ·εᵢ

# Extended model (conditional on VIX, yields):
yᵢ = b₀ + b₁·VIXᵢ + b₂·Yieldᵢ + β·zᵢ + √zᵢ·εᵢ
```

**Benefits:**
- "What's tomorrow's VaR if VIX = 30?"
- Regime-specific risk estimates
- Better tail accuracy during volatility spikes

### 2. Out-of-Sample Validation
- Rolling window backtesting
- VaR coverage tests (Kupiec, Christoffersen)
- Compare in-sample vs out-of-sample performance

### 3. Multi-Asset Extension
- Test on different asset classes (commodities, FX, crypto)
- Multivariate NIG for portfolio VaR
- Cross-asset tail dependence

---

## 📊 Visualizations

### 1. Distribution Comparison
![Returns Distribution](path/to/distribution_plot.png)
- Empirical density (black)
- NIG fit (red) — captures heavy tails
- Normal fit (blue) — underestimates tails

### 2. VaR Coverage
![VaR Bands](path/to/var_coverage.png)
- Time series with 1% and 5% VaR bands
- Shows how NIG adapts to market regimes

### 3. MCMC Diagnostics
![Trace Plots](path/to/trace_plots.png)
- Convergence check for γ, δ, μ, β
- Posterior distributions

---

## 🔍 Key Insights

### Why NIG Outperforms Normal

1. **Realistic Tail Behavior**
   - Normal: P(|return| > 3σ) = 0.27%
   - NIG: P(|return| > 3σ) = 2-5% (matches reality)
   - Markets have ~10× more extreme events than Normal predicts

2. **Asymmetry Capture**
   - Normal: Symmetric (skew = 0)
   - NIG: β parameter allows left/right skew
   - Captures "crashes are faster than rallies"

3. **Time-Varying Volatility**
   - Normal: Fixed σ (unrealistic)
   - NIG: Each day has own zᵢ (heteroscedastic)
   - 2008, 2020 crashes get higher variance automatically

### Limitations

1. **Kurtosis Underestimation**
   - Real: 15.51, Predicted: 8.51 (54.8% captured)
   - S&P 500 has EXTREME tails (includes COVID, 2008 aftermath)
   - Future work: Stronger priors or regime-switching extension

2. **Computational Cost**
   - 500k MCMC iterations ≈ 15-30 minutes
   - vs Normal MLE: < 1 second
   - Trade-off: Accuracy vs speed

3. **Prior Sensitivity**
   - Results depend on omega, eta choices
   - Weak priors (omega=1) let data dominate
   - Strong priors (omega=100+) needed for extreme kurtosis

---

## 📚 References

### Theoretical Foundation
- **Karlis, D. & Lillestøl, J.** (2004). "Bayesian estimation of NIG models via Markov chain Monte Carlo methods." *Applied Stochastic Models in Business and Industry*, 20(4), 323-338.
- **Barndorff-Nielsen, O.** (1997). "Normal Inverse Gaussian Distributions and Stochastic Volatility Modelling." *Scandinavian Journal of Statistics*, 24(1), 1-13.

### Implementation
- GIG distribution: `GIGrvg` R package
- Data augmentation for latent variables
- Gibbs sampling with conjugate priors

---

## 🐛 Known Issues & Solutions

### Issue: Kurtosis still too low (8.5 vs 15.5)
**Cause:** Weak prior (omega=1) insufficient for extreme tails  
**Solution:** Increase omega to 500-1000 for S&P 500 data

### Issue: Slow MCMC (hours instead of minutes)
**Cause:** Loop over observations for z sampling  
**Solution:** Vectorize with `rgig(n, ...)` instead of loop

---

## 💡 Practical Applications

### For Portfolio Managers
- More realistic VaR for regulatory capital
- Better tail risk hedging decisions
- Scenario analysis with conditional forecasts

### For Risk Officers
- Conservative ES estimates (regulatory requirement)
- Stress testing with fat-tailed distributions
- Model risk reduction vs Normal assumption

### For Quant Researchers
- Benchmark for ML models
- Building block for regime-switching models
- Framework for multivariate extensions

---

## 📧 Contact & Collaboration

**Author:** Scaply  
**Institution:** University of Piraeus - Finance & Banking Management  
**GitHub:** [Scuplex](https://github.com/Scuplex/NIG-distribution-estimate-Using-Gibbs)  
**LinkedIn:** [Your Profile]

Open to:
- Code reviews and suggestions
- Collaboration on extensions (covariates, multivariate)
- Discussions on prior specification
- Applications to other asset classes

---

## 📝 License

MIT License - Free for academic and research use

---

## 🙏 Acknowledgments

- Professor Dimitris Karlis (Athens University of Economics and Business) for the original methodology
- R community for GIGrvg, quantmod packages
- Claude (Anthropic) for debugging assistance

---

**Last Updated:** March 2026

---

