# Standardized Sudent-t distribution (mean 0, variance 1)
dts <- function(x, df, m = 0, s = 1) 
  dt((x - m) / (s * sqrt((df-2)/df)), df = df) / (s * sqrt((df-2)/df))
pts <- function(q, df, m = 0, s = s) pt((q - m) / (s * sqrt((df-2)/df)), df)
qts <- function(p, df) qt(p, df) / sqrt(df/(df-2))
rts <- function(n, df) rt(n, df) / sqrt(df/(df-2))

dgp_fun <- function(design, n, alpha, seed = NULL) {
  # Fix the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Get the parameters of the process
  if (design == 1) {
    x <- cbind(1, rchisq(n, 1))
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    gamma <- c(0, -1)
    eta <- c(1, 0)
  } else if (design == 2) {
    x <- cbind(1, rchisq(n, 1))
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    gamma <- c(0, -1)
    eta <- c(1, 0.5)
  } else if (design == 3) {
    cor <- 0.5
    rho <- 2 * sin(cor * pi/6)      # Pearson correlation
    P <- toeplitz(c(1, rho))        # Correlation matrix
    d <- nrow(P)                    # Dimension
    x <- cbind(1, pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P)))
    v <- rts(n, df = 5)
    q <- qts(alpha, df = 5)
    e <- integrate(qts, lower = 0, upper = alpha, df = 5)$value / alpha
    df <- function(x, m = 0, s = 1, df = 5) dt((x - m) / (s * sqrt((df-2)/df)), df = df) / (s * sqrt((df-2)/df))
    pf <- function(q, m = 0, s = 1, df = 5) pt((q - m) / (s * sqrt((df-2)/df)), df)
    closed_form <- 'student_t'
    gamma <- c(0, 1, -1)
    eta <- c(1, 1, 1)
  } else if (design == 4) {
    cor <- 0.5
    rho <- 2 * sin(cor * pi/6)      # Pearson correlation
    P <- toeplitz(c(1, rho))        # Correlation matrix
    d <- nrow(P)                    # Dimension
    x <- cbind(1, pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P)))
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    gamma <- c(0, -q, -2*e)
    eta <- c(1, 1, 2)
  }
  
  # Compute the dependent variable
  y <- x %*% gamma + (x %*% eta) * v
  
  # True parameters
  par_q <- gamma + q * eta
  par_e <- gamma + e * eta
  
  # Return results
  list(y = y, x = x, gamma = gamma, eta = eta, par_q = par_q, par_e = par_e,
       df = df, pf = pf, closed_form = closed_form)
}
