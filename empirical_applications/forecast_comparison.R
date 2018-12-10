rm(list = ls())
library(pacman)
p_load(xts, esreg, rugarch, RcppRoll, fGarch, reshape2, doParallel, dplyr, tidyr)
source('empirical_applications/forecast_comparison_functions.R')


# Load data ---------------------------------------------------------------

df <- read.csv("empirical_applications/data/realized_quantities/IBM")
date <- as.Date(as.character(df$X))

y <- xts(df$r_cc, order.by = date)  # oc vs cc?
x <- xts(df$rvol, order.by = date)

idx <- is.na(y) | is.na(x)
y <- y[!idx]
x <- x[!idx]

alpha <- 0.025


# Forecasting ---------------------------------------------------------------------------------

out <- do_predictions(y, x)


# Evaluation --------------------------------------------------------------

g1 <- 1
scores <- out %>% 
  group_by(model) %>% 
  spread(variable, value) %>% 
  summarise(
    hr = mean(r <= q),
    l1 = esr_loss(r=r, q=q, e=e, alpha=alpha, g1=g1, g2 = 1),
    l2 = esr_loss(r=r, q=q, e=e, alpha=alpha, g1=g1, g2 = 2),
    l3 = esr_loss(r=r, q=q, e=e, alpha=alpha, g1=g1, g2 = 3),
    l4 = esr_loss(r=r, q=q, e=e, alpha=alpha, g1=g1, g2 = 4),
    l5 = esr_loss(r=r, q=q, e=e, alpha=alpha, g1=g1, g2 = 5)
  )

sum_stat <- out %>% 
  group_by(model) %>% 
  spread(variable, value) %>% 
  summarise(
    mean_q = mean(q),
    var_q = var(q) * 1e5,
    mean_e = mean(e),
    var_e = var(e) * 1e5
  )


# murphy <- function(r, q, e, alpha, v, return_mean=TRUE) {
#   x <- (v <= e) * (1/alpha * (r <= q) * (q - r) - (q - v)) + (v <= r) * (r - v)
#   if (return_mean) {
#     mean(x)
#   } else {
#     x
#   }
# }
# murpy_diff <- function(r, q1, e1, q2, e2, alpha, v) {
#   m1 <- murphy(r=r, q=q1, e=e1, alpha=alpha, v=v, return_mean = FALSE)
#   m2 <- murphy(r=r, q=q2, e=e2, alpha=alpha, v=v, return_mean = FALSE)
#   m <- m1 - m2
#   c(mean(m), sd(m) / sqrt(length(m)))
# }
# 
# murphy(r=df$r, q=df$q, e=df$e, alpha=alpha, v=v)
# murpy_diff(r=df0$r, q1=df0$q, e1=df0$e, q2=df1$q, e2=df1$e, alpha=alpha, v=v)
# 
# 
# x <- seq(-0.2, 0.2, length.out = 1000)
# loss <- sapply(list(df0, df1, df2, df3), function(df) {
#   sapply(x, function(v) murphy(r=df$r, q=df$q, e=df$e, alpha=alpha, v=v))
# })
# 
# colnames(loss) <- names
# loss <- data.frame(x=x, loss)
# 
# write.csv(loss, file = "../../Plots/Empirical_Application/input/losses", row.names = FALSE)
# 
# 
# loss_diff <- lapply(list(df1, df2, df3), function(df) {
#   t(sapply(x, function(v) murpy_diff(r=df0$r, q1=df0$q, e1=df$e, q2=df$q, e2=df$e, alpha=alpha, v=v)))
# })
# mean_loss_diff <- do.call(cbind, lapply(loss_diff, "[",, 1))
# sd_loss_diff <- do.call(cbind, lapply(loss_diff, "[",, 2))
# colnames(mean_loss_diff) <- colnames(sd_loss_diff) <- names[-1]
# mean_loss_diff <- data.frame(x=x, mean_loss_diff)
# sd_loss_diff <- data.frame(x=x, sd_loss_diff)
# write.csv(mean_loss_diff, file = "../../Plots/Empirical_Application/input/mean_loss_diff", row.names = FALSE)
# write.csv(sd_loss_diff, file = "../../Plots/Empirical_Application/input/sd_loss_diff", row.names = FALSE)
# 
