rm(list = ls())
library(pacman)
p_load(esreg, esback, quantreg, xts, rugarch, quantmod)
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

split <- 1000
df_is <- df[1:split,]
df_os <- df[(split+1):nrow(df),]

fit <- esreg(r ~ q1 + q2 | e1 + e2, data = df_is, alpha = alpha)

bq <- coef(fit)[1:3]
be <- coef(fit)[4:6]
pq <- cbind(1, as.matrix(df_os)[,2:3]) %*% bq
pe <- cbind(1, as.matrix(df_os)[,4:5]) %*% be

p_new <- cbind(df_os, qc=pq, ec=pe)

scores <- expand.grid(
  g1 = c(1, 2),
  g2 = c(1, 2, 3, 4, 5),
  l1= NA,
  l2= NA,
  lc= NA
)

for (idx_row in seq_len(nrow(scores))) {
  g1 <- scores$g1[idx_row]
  g2 <- scores$g2[idx_row]
  
  scores$l1[idx_row] <- esr_loss(r = p_new$r, q = p_new$q1, e = p_new$e1, alpha = alpha, g1 = g1, g2 = g2)
  scores$l2[idx_row] <- esr_loss(r = p_new$r, q = p_new$q2, e = p_new$e2, alpha = alpha, g1 = g1, g2 = g2)
  scores$lc[idx_row] <- esr_loss(r = p_new$r, q = p_new$qc, e = p_new$ec, alpha = alpha, g1 = g1, g2 = g2)
}

l1 <- esr_loss(r = p_new$r, q = p_new$q1, e = p_new$e1, alpha = alpha)
l2 <- esr_loss(r = p_new$r, q = p_new$q2, e = p_new$e2, alpha = alpha)
lc <- esr_loss(r = p_new$r, q = p_new$qc, e = p_new$ec, alpha = alpha)

c(l1, l2, lc)

