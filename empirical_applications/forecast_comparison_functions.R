do_predictions <- function(y, x, win = 1000) {
  file <- 'empirical_applications/data/forecast_comparison.rds'
  if (file.exists(file)) {
    out <- readRDS(file)
  } else {
    #registerDoParallel(cores=1)
    registerDoSEQ()
    
    N <- length(y)
    horizon <- (win+1):N
 
    out <- foreach(t=(win+1):N, .errorhandling = 'pass') %dopar% {
      print(paste0(t, ' / ', N))
 
      # Store data of the current step
      y_ <- as.numeric(y[(t-win):(t-1)])
      x_ <- as.matrix(x[(t-win):(t-1),])
      yf <- as.numeric(y[t])
      d_ <- index(y)[t]
      
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
      
      # Model
      c('date', 'model', 'variable', 'value')
      
      c(d_, 'esr', 'q')
      
      d_
      'esr'
      'variable'
      
      # Return
      data.frame(
        date = rep(d_, 4*3),
        model = rep(c("esr", "garch", "har", "hs"), each=3),
        variable = rep(c('r', 'q', 'e'), 4),
        value = c(yf, q_esr, e_esr, 
                  yf, q_garch, e_garch, 
                  yf, q_har, e_har, 
                  yf, q_hs, e_hs),
        stringsAsFactors = FALSE
      )
    }
    out <- out[!sapply(out, function(x) inherits(x, "simpleError"))]
    out <- do.call('rbind', out)
    saveRDS(out, file)
  }
}

