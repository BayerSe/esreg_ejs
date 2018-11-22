library(esreg)
library(quantreg)
library(xts)
library(rugarch)


# Functions ---------------------------------------------------------------


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


# Main Code ---------------------------------------------------------------

p <- quantmod::getSymbols("^GSPC", auto.assign = FALSE, 
                          from = "2000-01-01", to = "2018-10-31")
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
  q[t] <- quantile(r_tmp, probs = alpha)
  e[t] <- mean(r_tmp[r_tmp <= q[t]])
  s[t] <- sd(r_tmp)
}
df2 <- data.frame(r = r, mu=NA, s = s, q = q, e = e)
df2 <- df2[-(1:(n0+250)),]


df <- data.frame(
  r = df1$GSPC.Adjusted,
  q1 = df1$q,
  e1 = df1$e,
  q2 = df2$q,
  e2 = df2$e
)

last_n <- 20
win <- 1000
index <- (nrow(df)-last_n+1):nrow(df)

pred <- lapply(index, function(idx) {
  df_ <- df[(idx-win):(idx-1),]
  fit <- esreg(r ~ q1 + q2 | e1 + e2, data = df_, alpha = alpha)  
  pq <- as.numeric(c(1, as.numeric(df[idx, c('q1', 'q2')])) %*% fit$coefficients_q)
  pe <- as.numeric(c(1, as.numeric(df[idx, c('e1', 'e2')])) %*% fit$coefficients_e)
  list(c(pq, pe), fit$coefficients)
})
p <- t(sapply(pred, function(x) x[[1]]))
coef <- t(sapply(pred, function(x) x[[2]]))

plot.ts(coef[,c(1,2,3)], plot.type = 'single')
plot.ts(coef[,c(4,5,6)], plot.type = 'single')


p_new <- cbind(df[index,], p)
plot.ts(p_new, plot.type = 'single')

esreg::esr_loss(r = p_new$r, q = p_new$q1, e = p_new$e1, alpha = alpha)
esreg::esr_loss(r = p_new$r, q = p_new$q2, e = p_new$e2, alpha = alpha)
esreg::esr_loss(r = p_new$r, q = p_new$`1`, e = p_new$`2`, alpha = alpha)
