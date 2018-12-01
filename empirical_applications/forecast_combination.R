library(esreg)
library(quantreg)
library(xts)
library(rugarch)
library(quantmod)


# Functions ---------------------------------------------------------------

#' @title Weighted Historical Simulation
#' @description Computes the Value-at-Risk and the Expected Shortfall
#' using the weighted historial simulation technique
#' @param x Return series, ordered by time
#' @param alpha Quantile of interest in [0, 1]
#' @param nu Geometric ecay factor in [0, 1]. 1 corresponds to standard
#' historial simulation
#' @return A vector with the predicted the Value-at-Risk and the Expected Shortfall
weighted_hs <- function(x, alpha, nu=1) {
  # Compute the weights
  n <- length(x)
  if (nu < 1) {
    w <- nu^(n:1-1) * (1-nu) / (1-nu^n)  
  } else {
    w <- rep(1/n, n)
  }
  
  # Empirical distribution
  ordered_x <- order(x)
  ecdf <- cbind(x[ordered_x], cumsum(w[ordered_x])) 
  
  # Compute the quantile and the expected shortfall
  s <- min(which(ecdf[,2] >= alpha))    
  q <- ecdf[s, 1]
  e <- mean(ecdf[1:s,1])
  s <- var(x)
  
  # Return estimates
  c(q, e, s)
}

filter_spec <- function(spec, r, burn_in, alpha) {
  
  # Filter the data with the provided specification
  filter <- ugarchfilter(spec, r)
  mu <- as.numeric(fitted(filter))
  s <- as.numeric(sigma(filter))
  
  # Quantile function of the innovations
  qf <- function(x) qdist(distribution=spec@model$modeldesc$distribution, p=x, mu=0, sigma=1,
                          shape=spec@model$fixed.pars$shape, skew=spec@model$fixed.pars$skew, 
                          lambda=spec@model$fixed.pars$ghlambda)
  cdf <- function(x) pdist(distribution=spec@model$modeldesc$distribution, q=x, mu=0, sigma=1,
                           shape=spec@model$fixed.pars$shape, skew=spec@model$fixed.pars$skew)
  
  # VaR and ES of the innovations
  vq <- qf(alpha)
  ve <- integrate(qf, 0, alpha)$value / alpha
  
  # VaR and ES of the returns
  q <- mu + s * vq
  e <- mu + s * ve
  
  # Returns filtered with the cdf
  x <- cdf((r - mu) / s)
  
  # Cummulative hit series
  H <- 1/alpha * (alpha - x) * (x <= alpha)   
  
  data.frame(r = r, mu = mu, s = s, q = q, e = e, H = H)[-(1:burn_in),]
}


# Main part ---------------------------------------------------------------

p <- quantmod::getSymbols("^GSPC", auto.assign = FALSE, 
                          from = "2000-01-01", to = "2018-11-30")
p <- p[,6]
r <- xts::diff.xts(p, log=TRUE, na.pad=FALSE) * 100

alpha <- 0.025
n0 <- 250


# Alternative 1: RiskMetrics
spec <- ugarchspec(mean.model=list(armaOrder=c(0, 0), include.mean=FALSE),
                   variance.model=list(model="iGARCH", garchOrder=c(1, 1)),
                   distribution.model="norm",
                   fixed.pars=list(omega=0, alpha1=1-0.94))
df1 <- filter_spec(spec = spec, r = r, burn_in = 250, alpha=alpha)
df1 <- df1[-(1:n0),]

# Alternative 2: Historical Simulation
q <- e <- s <- rep(NA, length(r))
for (t in (250+1):length(r)) {
  r_tmp <- r[(t-250):(t-1)]
  r_ <- weighted_hs(as.numeric(r_tmp), alpha = alpha, nu = 0.99)
  q[t] <- r_[1]
  e[t] <- r_[2]
  s[t] <- r_[3]
  
  #q[t] <- quantile(r_tmp, probs = alpha)
  #e[t] <- mean(r_tmp[r_tmp <= q[t]])
  #s[t] <- sd(r_tmp)
}
df2 <- data.frame(r = r, mu=NA, s = s, q = q, e = e)
df2 <- df2[-(1:(n0+250)),]


df <- data.frame(
  r = df1$GSPC.Adjusted,
  q1 = df1$q, e1 = df1$e,
  q2 = df2$q, e2 = df2$e
)

# In-sample combination ---------------------------------------------------

fit <- esreg(r ~ q1 + q2 | e1 + e2, data = df, alpha = alpha)
fit_r <- esreg(r ~ 1, data = df, alpha = alpha)

summary(fit, cond_var='scl_sp', sparsity='nid')



pq <- fitted(fit)[,1]
pe <- fitted(fit)[,2]
pqr <- fitted(fit_r)[,1]
per <- fitted(fit_r)[,2]

p_new <- cbind(df, qc=pq, ec=pe, qr=pqr, er=per)

l1 <- esr_loss(r = p_new$r, q = p_new$q1, e = p_new$e1, alpha = alpha)
l2 <- esr_loss(r = p_new$r, q = p_new$q2, e = p_new$e2, alpha = alpha)
lc <- esr_loss(r = p_new$r, q = p_new$qc, e = p_new$ec, alpha = alpha)
lr <- esr_loss(r = p_new$r, q = p_new$qr, e = p_new$er, alpha = alpha)

1 - c(l1, l2, lc) / lr

plot.ts(cbind(p_new$r, p_new$q1, p_new$e1), plot.type = 'single', col=1:3)
plot.ts(cbind(p_new$r, p_new$q2, p_new$e2), plot.type = 'single', col=1:3)
plot.ts(cbind(p_new$r, p_new$qc, p_new$ec), plot.type = 'single', col=1:3)

        