rm(list = ls())

library(esreg)
library(quantreg)
library(xts)

df <- read.csv('/home/sebastian/Seafile/Repos/github/RealizedQuantities/in/EURUSD_returns_1min.csv',
               row.names = NULL, header = FALSE) 
r <- apply(df, 2, sum)
rv <- apply(df, 2, function(x) sum(x^2))
rv_p <- apply(df, 2, function(x) sum(x^2 * (x > 0)))
rv_m <- apply(df, 2, function(x) sum(x^2 * (x < 0)))

rvol <- sqrt(rv)
rvol_p <- sqrt(rv_p)
rvol_m <- sqrt(rv_m)

df <- data.frame(r=r, rvol=rvol, rvol_p=rvol_p, rvol_m=rvol_m)

alpha <- 0.05

rq(r ~ rvol, tau = alpha, data = df)
rq(r ~ rvol_m + rvol_p, tau = alpha, data = df)

esreg(r ~ rvol, data = df, alpha = alpha)
esreg(r ~ rvol_m + rvol_p, data = df, alpha = alpha)


esreg(r ~ rvol, data = df, alpha = alpha)

