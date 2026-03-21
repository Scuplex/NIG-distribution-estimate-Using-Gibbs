# Data pull for the tickers 
library(writexl)
library(quantmod)

start_date <- "2007-05-01" 
end_date <- "2026-01-01"
tickers    <- c("SPY", "BIL", "VEU", "TLT") # ticker lsit 

getSymbols(tickers, from = start_date, to = end_date, periodicity = "daily")
prices <- merge(Ad(SPY), Ad(BIL), Ad(VEU), Ad(TLT)) # merge into one matrix
colnames(prices) <- tickers # give names to the columns of matrix

adj_close <- prices[endpoints(prices, "months")] # get end of month prices from the daily freq
adj_close <- na.omit(adj_close) # cause BIL didnt trade for most of May 2007

#save ino excel
adj_close_df <- data.frame(Date = index(adj_close), coredata(adj_close))
write_xlsx(adj_close_df, "prices_eom.xlsx")
