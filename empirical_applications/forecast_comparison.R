m(list = ls())
library(pacman)
p_load(xts, esreg, rugarch, RcppRoll, doParallel, fGarch, reshape2)


df <- read.csv("data/realized_quantities/IBM")
date <- as.Date(as.character(df$X))

y <- xts(df$r_cc, order.by = date)
x <- xts(df$rvol, order.by = date)

idx <- is.na(y) | is.na(x)
y <- y[!idx]
x <- x[!idx]

alpha <- 0.025



# Forecasting ---------------------------------------------------------------------------------

registerDoParallel(cores=4)

win <- 1000
N <- length(y)
horizon <- (win+1):N
plot(y[horizon])


out <- foreach (t=(win+1):N) %dopar% {
  
  # Store data of the current step
  y_ <- as.numeric(y[(t-win):(t-1)])
  x_ <- as.matrix(x[(t-win):(t-1),])
  yf <- as.numeric(y[t])
  
  # Expected Shortfall Regression
  fit_esr <- esreg(y_[2:win] ~ x_[1:(win-1)], alpha=alpha)
  q_esr <- as.numeric(c(1, x_[win,]) %*% fit_esr$coefficients_q)
  e_esr <- as.numeric(c(1, x_[win,]) %*% fit_esr$coefficients_e)
  
  # GARCH
  fit_garch <- garchFit(~arma(1,0) + garch(1,1), data=y_, cond.dist="std", trace=FALSE)
  f <- function(x, fit) qstd(p=x, mean=0, sd=1, nu=fit@fit$par["shape"]) 
  fcst_garch <- predict(fit_garch, n.ahead=1)
  q_garch <- fcst_garch$meanForecast + fcst_garch$standardDeviation * f(alpha, fit_garch)
  e_garch <- fcst_garch$meanForecast + fcst_garch$standardDeviation * integrate(f, 0, alpha, fit=fit_garch)$value / alpha
  
  # HAR
  data_har <- cbind(x_, 
                    roll_meanr(x_, n=5, fill=NA), 
                    roll_meanr(x_, n=22, fill=NA))
  fit_har <- lm(x_[2:win] ~ data_har[1:(win-1),])
  
  s_f <- fit_har$coefficients %*% c(1, tail(data_har, 1))
  q_har <- s_f * qnorm(alpha)
  e_har <- s_f * -dnorm(qnorm(alpha)) / alpha
  
  # Historical simulation
  q_hs <- quantile(tail(y_, 250), alpha)
  e_hs <- mean(tail(y_, 250)[tail(y_, 250) <= q_hs])
  
  # Return everything as a list
  list(esr     = unname(c(yf, q_esr, e_esr)),
       garch   = unname(c(yf, q_garch, e_garch)),
       har     = unname(c(yf, q_har, e_har)),
       hs      = unname(c(yf, q_hs, e_hs))
  )
}

names <- c("esr", "garch", "har", "hs")
df0 <- as.data.frame(t(sapply(out, "[[", "esr")))
df1 <- as.data.frame(t(sapply(out, "[[", "garch")))
df2 <- as.data.frame(t(sapply(out, "[[", "har")))
df3 <- as.data.frame(t(sapply(out, "[[", "hs")))
names(df0) <- names(df1) <- names(df2) <- names(df3) <- c("r", "q", "e")

q <- cbind(df0$q, df1$q, df2$q, df3$q)
e <- cbind(df0$e, df1$e, df2$e, df3$e)

apply(q, 2, mean)
apply(q, 2, var) * 1e5

apply(e, 2, mean)
apply(e, 2, var) * 1e5

plot.zoo(q, plot.type = "single", col=1:3)
plot.zoo(e, plot.type = "single", col=1:3)

loss <- sapply(list(df0, df1, df2, df3), function(df) {
  c(
    mean(df$r <= df$q),
    gpl(r=df$r, q=df$q, alpha=alpha) * 1e4,
    esr_loss(r=df$r, q=df$q, e=df$e, alpha=alpha, g2 = 1),
    esr_loss(r=df$r, q=df$q, e=df$e, alpha=alpha, g2 = 2),
    esr_loss(r=df$r, q=df$q, e=df$e, alpha=alpha, g2 = 3),
    esr_loss(r=df$r, q=df$q, e=df$e, alpha=alpha, g2 = 4),
    esr_loss(r=df$r, q=df$q, e=df$e, alpha=alpha, g2 = 5)
  )
})
colnames(loss) <- names
rownames(loss) <- c("hr", "tl", "jl1", "jl2", "jl3", "jl4", "jl5")
loss

esreg(df0$r ~ df0$e, alpha=alpha)
esreg(df1$r ~ df1$e, alpha=alpha)
esreg(df2$r ~ df2$e, alpha=alpha)


murphy <- function(r, q, e, alpha, v, return_mean=TRUE) {
  x <- (v <= e) * (1/alpha * (r <= q) * (q - r) - (q - v)) + (v <= r) * (r - v)
  if (return_mean) {
    mean(x)
  } else {
    x
  }
}
murpy_diff <- function(r, q1, e1, q2, e2, alpha, v) {
  m1 <- murphy(r=r, q=q1, e=e1, alpha=alpha, v=v, return_mean = FALSE)
  m2 <- murphy(r=r, q=q2, e=e2, alpha=alpha, v=v, return_mean = FALSE)
  m <- m1 - m2
  c(mean(m), sd(m) / sqrt(length(m)))
}

murphy(r=df$r, q=df$q, e=df$e, alpha=alpha, v=v)
murpy_diff(r=df0$r, q1=df0$q, e1=df0$e, q2=df1$q, e2=df1$e, alpha=alpha, v=v)


x <- seq(-0.2, 0.2, length.out = 1000)
loss <- sapply(list(df0, df1, df2, df3), function(df) {
  sapply(x, function(v) murphy(r=df$r, q=df$q, e=df$e, alpha=alpha, v=v))
})

colnames(loss) <- names
loss <- data.frame(x=x, loss)

write.csv(loss, file = "../../Plots/Empirical_Application/input/losses", row.names = FALSE)


loss_diff <- lapply(list(df1, df2, df3), function(df) {
  t(sapply(x, function(v) murpy_diff(r=df0$r, q1=df0$q, e1=df$e, q2=df$q, e2=df$e, alpha=alpha, v=v)))
})
mean_loss_diff <- do.call(cbind, lapply(loss_diff, "[",, 1))
sd_loss_diff <- do.call(cbind, lapply(loss_diff, "[",, 2))
colnames(mean_loss_diff) <- colnames(sd_loss_diff) <- names[-1]
mean_loss_diff <- data.frame(x=x, mean_loss_diff)
sd_loss_diff <- data.frame(x=x, sd_loss_diff)
write.csv(mean_loss_diff, file = "../../Plots/Empirical_Application/input/mean_loss_diff", row.names = FALSE)
write.csv(sd_loss_diff, file = "../../Plots/Empirical_Application/input/sd_loss_diff", row.names = FALSE)

