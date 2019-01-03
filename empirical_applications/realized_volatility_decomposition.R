rm(list = ls())

library(esreg)
library(quantreg)
library(xts)


df <- read.csv('https://github.com/BayerSe/RealizedQuantities/raw/master/out/realized_quantities_IBM_cts.csv')

subs <- as.Date(df$X) > as.Date('2004-01-01')
df <- df[subs,]

df_ <- cbind(df$rvol_p, df$rvol_m)
colMeans(df_) * 10^5
cor(df_)

alpha <- 0.10

rq(r_oc ~ rvol, tau = alpha, data = df)
rq(r_oc ~ rvol_m + rvol_p, tau = alpha, data = df)
rq(r_oc ~ ivol + jvol, tau = alpha, data = df)

esreg(r_oc ~ rvol, data = df, alpha = alpha)
esreg(r_oc ~ rvol_m + rvol_p, data = df, alpha = alpha)
esreg(r_oc ~ ivol + jvol, data = df, alpha = alpha)

