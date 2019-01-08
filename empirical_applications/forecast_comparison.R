rm(list = ls())
library(pacman)
p_load(xts, esreg, fGarch, doParallel, dplyr, tidyr, xtable, RcppRoll, reticulate)
use_virtualenv('/home/sebastian/.virtualenvs/esreg_theory/') 
source_python('../esreg_python_functions/main.py')
source('empirical_applications/empirical_application_functions.R')

tbl_folder <- '../../Revision EJS/Paper/v1/tables/'
img_folder <- '../../Revision EJS/Paper/v1/plots/'


# Functions ---------------------------------------------------------------

do_predictions <- function(y, x, alpha, win = 1000) {
  file <- 'empirical_applications/data/forecast_comparison.rds'
  if (file.exists(file)) {
    out <- readRDS(file)
  } else {
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
      
      # Return
      data_frame(
        date = rep(d_, 4*3),
        model = rep(c("esr", "garch", "har", "hs"), each=3),
        variable = rep(c('r', 'q', 'e'), 4),
        value = c(yf, q_esr, e_esr, 
                  yf, q_garch, e_garch, 
                  yf, q_har, e_har, 
                  yf, q_hs, e_hs)
      )
    }
    out <- out[!sapply(out, function(x) inherits(x, "simpleError"))]
    out <- do.call('rbind', out)
    saveRDS(out, file)
  }
}


# Load data ---------------------------------------------------------------

df <- read.csv("empirical_applications/data/realized_quantities/IBM")
date <- as.Date(as.character(df$X))

y <- xts(df$r_cc, order.by = date) * 100
x <- xts(df$rvol, order.by = date) * 100

idx <- is.na(y) | is.na(x)
y <- y[!idx]
x <- x[!idx]

alpha <- 0.025

out <- do_predictions(y = y, x = x, alpha = alpha, win = 1000)


# Evaluation --------------------------------------------------------------

df_os <- out %>%
  unite(temp, model, variable) %>%
  spread(temp, value) %>%
  rename(r = esr_r) %>%
  select(-contains('_r'))

scores <- expand.grid(
  g1 = c(1, 2),
  g2 = c(1, 2, 3, 4, 5),
  l_esr = NA, p_esr = NA,
  l_hs = NA, p_hs = NA,
  l_garch = NA, p_garch = NA,
  l_har = NA, p_har = NA
)

for (idx_row in seq_len(nrow(scores))) {
  print(idx_row / nrow(scores))
  g1 <- scores$g1[idx_row]
  g2 <- scores$g2[idx_row]
  
  l_esr <- esr_loss(r = df_os$r, q = df_os$esr_q, e = df_os$esr_e, 
                   alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  l_hs <- esr_loss(r = df_os$r, q = df_os$hs_q, e = df_os$hs_e, 
                   alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  l_har <- esr_loss(r = df_os$r, q = df_os$har_q, e = df_os$har_e, 
                     alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  l_garch <- esr_loss(r = df_os$r, q = df_os$garch_q, e = df_os$garch_e, 
                    alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)  
  loss <- data.frame(ESR=l_esr, HS=l_hs, HAR=l_har, GARCH=l_garch)
  
  # Compute scores
  scores$l_esr[idx_row] <- mean(l_esr)
  scores$l_hs[idx_row] <- mean(l_hs)
  scores$l_garch[idx_row] <- mean(l_garch)
  scores$l_har[idx_row] <- mean(l_har)
  
  # Compute MCS
  mcs <- compute_mcs(loss=loss, reps=10000L, bootstrap='stationary', 
                     statistic='R', block_size=10L)
  scores$p_esr[idx_row] <- mcs['ESR', 'Pvalue']
  scores$p_hs[idx_row] <- mcs['HS', 'Pvalue']
  scores$p_garch[idx_row] <- mcs['GARCH', 'Pvalue']
  scores$p_har[idx_row] <- mcs['HAR', 'Pvalue']
}

# Table -------------------------------------------------------------------


get_table <- function(g1_, conflevel=0.1) {
  p <- scores %>% subset(g1 == g1_) %>% select(c(p_esr, p_hs, p_garch, p_har)) %>% t()
  l <- scores %>% subset(g1 == g1_) %>% select(c(l_esr, l_hs, l_garch, l_har)) %>% t()
  
  col <- function(loss, pval, conflevel=0.1){
    ifelse(pval > conflevel, 
           ifelse(loss == min(loss), sprintf("%.3f^{**}", loss), sprintf("%.3f^{*}", loss)),
           sprintf("%.3f", loss))} 
  tab <- sapply(1:ncol(l), function(i) col(l[,i], p[,i]))
  tab <- cbind(c('ESR', 'HS', 'GARCH', 'HAR'), tab)
  g1_fun <- ifelse(g1_==1, '$G_1(z) = z$', '$G_1(z) = 0$')
  tab <- cbind(c(paste0('\\multirow{4}{*}{', g1_fun, '}'), '', '', ''), tab)
  
  tab
}

tab <- rbind(get_table(1), get_table(2))
print(xtable(tab), file = paste0(tbl_folder, 'forecast_comparison_scores.txt'),
      include.rownames = FALSE, include.colnames = FALSE, 
      sanitize.text.function = function(x) x, booktabs = TRUE,
      comment = FALSE, only.contents = TRUE, hline.after = c(4,8))


# Murphy diagram ----------------------------------------------------------


x <- seq(-10, 0, length.out = 1000)
models <- unique(out$model)[unique(out$model) != 'esr']

loss_diff <- lapply(models, function(est) {
  df1 <- out %>% subset(model == est) %>% spread(variable, value)
  df0 <- out %>% subset(model == 'esr') %>% spread(variable, value)
  
   t(sapply(x, function(v) murpy_diff(r=df0$r, q1=df0$q, e1=df0$e, 
                                      q2=df1$q, e2=df1$e, alpha=alpha, v=v)))
 })

mean_loss_diff <- data.frame(do.call(cbind, lapply(loss_diff, "[",, 1)))
sd_loss_diff <- data.frame(do.call(cbind, lapply(loss_diff, "[",, 2)))
colnames(mean_loss_diff) <- colnames(sd_loss_diff) <- toupper(models)
mean_loss_diff <- mean_loss_diff[,c('HS', 'GARCH', 'HAR')]
sd_loss_diff <- sd_loss_diff[,c('HS', 'GARCH', 'HAR')]

mean_loss_diff$x <- x
sd_loss_diff$x <- x

plot_murphy_diagram(loss_mean=mean_loss_diff, loss_sd=sd_loss_diff, cl=0.9, 
                    file=paste0(img_folder, 'forecast_comparison_murphy.pdf'))
