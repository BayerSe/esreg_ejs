library(esreg)
library(quantreg)

df <- read.csv('https://github.com/BayerSe/RealizedQuantities/raw/master/out/realized_quantities_IBM_cts.csv')

r <- df$r_oc[2:nrow(df)]
rvol <- df$rvol[2:nrow(df)]
rvol_m <- df$rvol_m[2:nrow(df)]
rvol_p <- df$rvol_p[2:nrow(df)]
print(cor(cbind(rvol, rvol_m, rvol_p)))

alpha <- 0.05

fit_rq <- rq(r ~ rvol, tau = alpha)
fit <- esreg(r ~ rvol, alpha = alpha)

fit_rq <- rq(r ~ rvol_m + rvol_p, tau = alpha)
fit <- esreg(r ~ rvol_m + rvol_p, alpha = alpha)
