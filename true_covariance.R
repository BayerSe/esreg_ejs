rm(list = ls())

library(esreg)
library(Rcpp)

source("dgp_functions.R")
sourceCpp("cpp_functions.cpp")

path <- "data/true_covariance/"
dir.create(path, showWarnings = FALSE, recursive = TRUE)

# Simulation setup
setup <- expand.grid(
  alpha   = 0.025,
  n       = 10^9,
  design  = c(1, 2, 3, 4),
  g1      = c(1, 2),
  g2      = c(1, 2, 3, 4, 5)
)

# Compute the covariance for all settings
for (idx in 1:nrow(setup)) {
  st <- setup[idx,]
  
  # File name
  file <- paste0(path, "design_", st$design, "_g1_", st$g1, "_g2_", st$g2, ".rds")
  if (file.exists(file)) next
  file.create(file)
  print(file)
  
  # Simulate data and extract some components
  df <- dgp_fun(design = st$design, n = st$n, alpha = st$alpha, seed = 1)
  
  # Scale data
  if (st$g2 %in% c(1, 2, 3)) {
    max_y <- max(df$y)
    df$y <- df$y - max_y
    df$par_q[1] <- df$par_q[1] - max_y
    df$par_e[1] <- df$par_e[1] - max_y
    df$gamma[1] <- df$gamma[1] - max_y
  }
  
  # True quantile, shortfall, mean and variance of Y and the error term
  xq <- as.numeric(df$x %*% df$par_q)
  xe <- as.numeric(df$x %*% df$par_e)
  mu_y <- as.numeric(df$x %*% df$gamma)
  sd_y <- as.numeric(df$x %*% df$eta)
  mu_u <- mu_y - xq
  sd_u <- sd_y
  
  # G-Function values
  G1_prime_xq <- G_vec(z = xq, g = "G1_prime", type = st$g1)
  G2_xe <- G_vec(z = xe, g = "G2", type = st$g2)
  G2_prime_xe <- G_vec(z = xe, g = "G2_prime", type = st$g2)
  
  # True density
  dens <- df$df(x = xq, m = mu_y, s = sd_y)
  
  # Truncated conditional variance
  if (df$closed_form == "norm") {
    b <- - mu_u / sd_u
    cond_var <- sd_u^2 * (1 - b * dnorm(b) / pnorm(b) - (dnorm(b) / pnorm(b))^2)
  } else if (df$closed_form == "student_t") {
    dof <- 5
    t1 <- (dof - 2) / dof
    beta <- - mu_u / sd_u
    if (dof < 300) {
      k <- gamma((dof+1)/2) / gamma(dof/2) / sqrt(dof*pi) / stats::pt(beta, df = dof)
    } else {
      k <- sqrt(dof/2) / sqrt(dof*pi) / stats::pt(beta, df = dof)
    }
    m1 <- k * dof / (dof - 1) * (-(1+beta^2/dof)^(-(dof-1)/2))
    m2 <- (dof - 1) / t1 * (stats::pt(beta*sqrt(t1), df = dof-2) / stats::pt(beta, df = dof)) - dof
    cond_var <- t1 * sd_u^2 * (m2 - m1^2)
  } else {
    cond_var <- rep(NA, n)
    for (i in 1:n) {
      # Conditional truncated distribution
      tdf <- function(x) df(x, m = mu_u[i], s = sd_u[i]) / pf(0, m = mu_u[i], s = sd_u[i])
      # Variance
      m1 <- integrate(function(x) x * tdf(x),   lower = -Inf, upper = 0)$value
      m2 <- integrate(function(x) x^2 * tdf(x), lower = -Inf, upper = 0)$value
      cond_var[i] <- m2 - m1^2
    }
  }
  
  # Compute the covariance
  cov <- esreg:::l_esreg_covariance(xq = df$x, xe = df$x, xbq = xq, xbe = xe, 
                                    G1_prime_xq = G1_prime_xq, G2_xe = G2_xe, 
                                    G2_prime_xe = G2_prime_xe, density = dens, 
                                    conditional_variance = cond_var, 
                                    alpha = st$alpha) * st$n
  
  cov_qr <- l_qr_covariance(x = df$x, density = dens, alpha = st$alpha) * st$n
  
  # Store as list
  out <- list(design = st$design, g1 = st$g1, g2 = st$g2, 
              cov = cov, cov_qr = cov_qr)
  print(out$cov)
  
  saveRDS(out, file = file)
}

