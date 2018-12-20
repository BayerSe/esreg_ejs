rm(list = ls())
library(pacman)
p_load(esreg, quantreg, xts, quantmod, xtable, dplyr, tidyr, reticulate)
use_virtualenv('/home/sebastian/.virtualenvs/esreg_theory/') 
source_python('../esreg_python_functions/main.py')
source('empirical_applications/empirical_application_functions.R')

tbl_folder <- '../../Revision EJS/Paper/v1/tables/'
img_folder <- '../../Revision EJS/Paper/v1/plots/'

# Functions ---------------------------------------------------------------

load_returns <- function(symbol='^GSPC', from='2000-01-01', to='2018-11-30') {
  p <- getSymbols(symbol, auto.assign = FALSE, from = from, to = to)
  p <- p[,6]
  r <- xts::diff.xts(p, log=TRUE, na.pad=FALSE) * 100
  unname(r)
}

mod_hs <- function(r, alpha, win = 250) {
  q <- e <- s <- rep(NA, length(r))
  for (t in (win+1):length(r)) {
    r_tmp <- r[(t-win):(t-1)]
    q[t] <- quantile(r_tmp, probs = alpha)
    e[t] <- mean(r_tmp[r_tmp <= q[t]])
  }
  
  df <- data_frame(
    date = rep(index(r), 3),
    model = 'HS',
    variable = rep(c('r', 'q', 'e'), each=length(r)),
    value = c(as.numeric(r), q, e)
  ) %>% arrange(date)
  df
}

mod_rm <- function(r, alpha, win=250) {
  a <- 0.96
  n <- length(r)
  s2 <- rep(NA, n)
  s2[1] <- var(r)
  for (i in 2:n) {
    s2[i] <- a * s2[i-1] + (1-a) * r[i-1]^2 
  }
  s <- sqrt(s2)
  q <- qnorm(alpha) * s
  e <- -dnorm(qnorm(alpha)) / alpha * s
  
  df <- data_frame(
    date = rep(index(r), 3),
    model = 'RM',
    variable = rep(c('r', 'q', 'e'), each=length(r)),
    value = c(as.numeric(r), q, e)
  ) %>% arrange(date)
  df
}


# Main part ---------------------------------------------------------------

alpha <- 0.025
r <- load_returns()
df_rm <- mod_rm(r, alpha)
df_hs <- mod_hs(r, alpha)
df <- rbind(df_rm, df_hs)

complete_obs <- df %>% 
  drop_na() %>% 
  subset(variable == 'q') %>% 
  group_by(date) %>% 
  count() %>% 
  subset(n == 2) %>% 
  pull(date)


# Forecast Combination ----------------------------------------------------

# Separate in-sample and out-of-sample data
split <- 1500

df_is <- df %>% 
  subset(date %in% complete_obs[1:split]) %>% 
  unite(temp, model, variable) %>% 
  spread(temp, value) %>% 
  rename(r = HS_r) 

df_os <- df %>% 
  subset(date %in% complete_obs[-c(1:split)]) %>% 
  unite(temp, model, variable) %>% 
  spread(temp, value) %>% 
  rename(r = HS_r) %>% 
  select(-contains('_r'))

# Do the actual combination
fit <- esreg(r ~ RM_q + HS_q | RM_e + HS_e, data = df_is, alpha = alpha)
pred_comb <- predict(fit, df_os)

# Append the combined forecasts
df_comb <- data_frame(
  date = rep(df_os$date, 3),
  model = 'Comb',
  variable = rep(c('r', 'q', 'e'), each=length(df_os$date)),
  value = c(df_os$r, pred_comb[,1], pred_comb[,2])
) %>% arrange(date)

df <- rbind(df, df_comb) %>% arrange(date)

# Update the out-of-sample forecasts
df_os <- df %>%
  subset(date %in% complete_obs[-c(1:split)]) %>%
  unite(temp, model, variable) %>%
  spread(temp, value) %>%
  rename(r = HS_r) %>%
  select(-contains('_r'))


# Evaluation --------------------------------------------------------------

scores <- expand.grid(
  g1 = c(1, 2),
  g2 = c(1, 2, 3, 4, 5),
  l_rm = NA, p_rm = NA,
  l_hs = NA, p_hs = NA,
  l_comb = NA, p_comb = NA
)

for (idx_row in seq_len(nrow(scores))) {
  print(idx_row / nrow(scores))
  g1 <- scores$g1[idx_row]
  g2 <- scores$g2[idx_row]
  
  l_rm <- esr_loss(r = df_os$r, q = df_os$RM_q, e = df_os$RM_e, 
                   alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  l_hs <- esr_loss(r = df_os$r, q = df_os$HS_q, e = df_os$HS_e, 
                   alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  l_comb <- esr_loss(r = df_os$r, q = df_os$Comb_q, e = df_os$Comb_e, 
                     alpha = alpha, g1 = g1, g2 = g2, return_mean = FALSE)
  loss <- data.frame(RiskMetrics=l_rm, Historical_Simulation=l_hs, 
                     Forecast_Combination=l_comb)
  
  # Compute scores
  scores$l_rm[idx_row] <- mean(l_rm)
  scores$l_hs[idx_row] <- mean(l_hs)
  scores$l_comb[idx_row] <- mean(l_comb)
  
  # Compute MCS
  mcs <- compute_mcs(loss=loss, reps=1000L, bootstrap='stationary', 
                     statistic='R', block_size=10L)
  scores$p_rm[idx_row] <- mcs['RiskMetrics', 'Pvalue']
  scores$p_hs[idx_row] <- mcs['Historical_Simulation', 'Pvalue']
  scores$p_comb[idx_row] <- mcs['Forecast_Combination', 'Pvalue']
}



# Table -------------------------------------------------------------------


get_table <- function(g1_) {
  p <- scores %>% filter(g1 == g1_) %>% select(c(p_rm, p_hs, p_comb)) %>% t()
  l <- scores %>% filter(g1 == g1_) %>% select(c(l_rm, l_hs, l_comb)) %>% t()
  
  col <- function(loss, pval, conflevel=0.1){
    ifelse(pval > conflevel, 
           ifelse(loss == min(loss), sprintf("%.3f^{**}", loss), sprintf("%.3f^{*}", loss)),
           sprintf("%.3f", loss))} 
  tab <- sapply(1:ncol(l), function(i) col(l[,i], p[,i]))
  tab <- cbind(c('RiskMetrics', 'Hist. Sim.', 'Combination'), tab)
  g1_fun <- ifelse(g1_==1, '$G_1(z) = z$', '$G_1(z) = 0$')
  tab <- cbind(c(paste0('\\multirow{3}{*}{', g1_fun, '}'), '', '', ''), tab)
  
  tab
}

tab <- rbind(get_table(1), get_table(2))
print(xtable(tab), file = paste0(tbl_folder, 'forecast_combination_scores.txt'),
      include.rownames = FALSE, include.colnames = FALSE, 
      sanitize.text.function = function(x) x, booktabs = TRUE,
      comment = FALSE, only.contents = TRUE, hline.after = c(3,6))


# Murphy diagram ----------------------------------------------------------


out <- df %>% subset(date %in% complete_obs[-c(1:(split+1))])
models <- unique(out$model)[unique(out$model) != 'Comb']

x <- seq(-10, 0, length.out = 1000)
loss_diff <- lapply(models, function(est) {
  df1 <- out %>% filter(model == est) %>% spread(variable, value)
  df0 <- out %>% filter(model == 'Comb') %>% spread(variable, value)
  
  t(sapply(x, function(v) murpy_diff(r=df0$r, q1=df0$q, e1=df0$e, 
                                     q2=df1$q, e2=df1$e, alpha=alpha, v=v)))
})

mean_loss_diff <- data.frame(do.call(cbind, lapply(loss_diff, "[",, 1)))
sd_loss_diff <- data.frame(do.call(cbind, lapply(loss_diff, "[",, 2)))
colnames(mean_loss_diff) <- colnames(sd_loss_diff) <- models
mean_loss_diff$x <- x
sd_loss_diff$x <- x

plot_murphy_diagram(loss_mean=mean_loss_diff, loss_sd=sd_loss_diff, cl=0.9, 
                    file=paste0(img_folder, 'forecast_combination_murphy.pdf'))
