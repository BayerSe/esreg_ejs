rm(list = ls())

library(esreg)
library(quantreg)
library(xts)
library(rugarch)
library(quantmod)

source('empirical_applications/forecast_combination_functions.R')

alpha <- 0.025
r <- load_returns()
df1 <- mod_riskmetrics(r, alpha)
df2 <- mod_weighted_hs(r, alpha)

df <- data.frame(
  r = df1$GSPC.Adjusted,
  q1 = df1$q, e1 = df1$e,
  q2 = df2$q, e2 = df2$e
)

# Forecast Combination ----------------------------------------------------

split <- 2000
df_is <- df[1:split,]
df_os <- df[(split+1):nrow(df),]

fit <- esreg(r ~ q1 + q2 | e1 + e2, data = df_is, alpha = alpha)
fit_r <- esreg(r ~ 1, data = df_is, alpha = alpha)

summary(fit, cond_var='scl_sp', sparsity='nid')

bqr <- coef(fit_r)
bq <- coef(fit)[1:3]
be <- coef(fit)[4:6]

pq <- cbind(1, as.matrix(df_os)[,2:3]) %*% bq
pe <- cbind(1, as.matrix(df_os)[,4:5]) %*% be
pqr <- rep(bqr[1], nrow(df_os))
per <- rep(bqr[2], nrow(df_os))

p_new <- cbind(df_os, qc=pq, ec=pe, qr=pqr, er=per)

l1 <- esr_loss(r = p_new$r, q = p_new$q1, e = p_new$e1, alpha = alpha)
l2 <- esr_loss(r = p_new$r, q = p_new$q2, e = p_new$e2, alpha = alpha)
lc <- esr_loss(r = p_new$r, q = p_new$qc, e = p_new$ec, alpha = alpha)
lr <- esr_loss(r = p_new$r, q = p_new$qr, e = p_new$er, alpha = alpha)

1 - c(l1, l2, lc) / lr

plot.ts(cbind(p_new$r, p_new$q1, p_new$e1), plot.type = 'single', col=1:3)
plot.ts(cbind(p_new$r, p_new$q2, p_new$e2), plot.type = 'single', col=1:3)
plot.ts(cbind(p_new$r, p_new$qc, p_new$ec), plot.type = 'single', col=1:3)

        