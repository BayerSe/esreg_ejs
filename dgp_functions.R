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
    z <- cbind(1, rchisq(n, 1))
    xq <- z
    xe <- z
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    gamma <- c(0, -1)
    eta <- c(1, 0)
    y <- z %*% gamma + (z %*% eta) * v  # Compute the dependent variable
    par_q <- gamma + q * eta  # True parameters
    par_e <- gamma + e * eta    
  } else if (design == 2) {
    z <- cbind(1, rchisq(n, 1))
    xq <- z
    xe <- z
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    gamma <- c(0, -1)
    eta <- c(1, 0.5)
    y <- z %*% gamma + (z %*% eta) * v  # Compute the dependent variable
    par_q <- gamma + q * eta  # True parameters
    par_e <- gamma + e * eta    
  } else if (design == 3) {
    cor <- 0.5
    rho <- 2 * sin(cor * pi/6)      # Pearson correlation
    P <- toeplitz(c(1, rho))        # Correlation matrix
    d <- nrow(P)                    # Dimension
    z <- cbind(1, pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P)))
    xq <- z
    xe <- z
    v <- rts(n, df = 5)
    q <- qts(alpha, df = 5)
    e <- integrate(qts, lower = 0, upper = alpha, df = 5)$value / alpha
    df <- function(x, m = 0, s = 1, df = 5) dt((x - m) / (s * sqrt((df-2)/df)), df = df) / (s * sqrt((df-2)/df))
    pf <- function(q, m = 0, s = 1, df = 5) pt((q - m) / (s * sqrt((df-2)/df)), df)
    closed_form <- 'student_t'
    gamma <- c(0, 1, -1)
    eta <- c(1, 1, 1)
    y <- z %*% gamma + (z %*% eta) * v  # Compute the dependent variable
    par_q <- gamma + q * eta  # True parameters
    par_e <- gamma + e * eta    
  } else if (design == 4) {
    cor <- 0.5
    rho <- 2 * sin(cor * pi/6)      # Pearson correlation
    P <- toeplitz(c(1, rho))        # Correlation matrix
    d <- nrow(P)                    # Dimension
    z <- cbind(1, pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P)))
    xq <- z[,c(1,3)]
    xe <- z[,c(1,2)]
    v <- rnorm(n)
    q <- qnorm(alpha)
    e <- -dnorm(qnorm(alpha)) / alpha
    df <- function(x, m = 0, s = 1) dnorm(x, mean = m, sd = s)
    pf <- function(q, m = 0, s = 1) pnorm(q, mean = m, sd = s)
    closed_form <- 'norm'
    #gamma <- c(0, -q, -2*e)
    #eta <- c(1, 1, 2)
    gamma <- c(0, -q, -e)
    eta <- c(0, 1, 1)
    y <- z %*% gamma + (z %*% eta) * v  # Compute the dependent variable
    par_q <- gamma[-2] + q * eta[-2]  # True parameters, and remove the par
    par_e <- gamma[-3] + e * eta[-3]  # which we do not use in the estimation
  }
  
  list(
    y = y, z = z, xq = xq, xe = xe,
    gamma = gamma, eta = eta, par_q = par_q, par_e = par_e,
    df = df, pf = pf, closed_form = closed_form
  )
}
