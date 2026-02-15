# Using linear regression
library(quantmod)
library(GIGrvg)
library(fBasics)

# Get SPY data for the MCMC
getSymbols("SPY", from = "2014-12-01", to = "2025-12-31", periodicity = "daily")
SPY_adj <- Ad(SPY)
R_t <- (SPY_adj / lag(SPY_adj)) - 1
y <- R_t[-1] # Remove NA