rm(list=ls())

library(pacman)
p_load(naturalsort, dplyr, reshape2, tidyr, xtable)


# Functions ---------------------------------------------------------------


fnorm <- function(mat, subset=NULL) {
  if (!is.null(subset)) {
    mat <- mat[subset, subset]  
  }   
  sqrt(mean((mat[lower.tri(mat, diag=TRUE)])^2))
}

collect_all_results <- function(sim_path) {
  dir_list <- list.files(paste0(sim_path, 'monte_carlo'), full.names = TRUE)
  dir_list <- naturalsort(dir_list)
  
  all_results <- list()
  for (idx in 1:length(dir_list)) {
    print(basename(dir_list[idx]))
    
    # Load results
    file_list <- list.files(dir_list[idx], full.names = TRUE)
    file_list <- file_list[file.size(file_list) > 0]
    
    #chunk <- as.numeric(sapply(file_list, 
    #                           function(s) gsub('.rds', '', strsplit(strsplit(s, '/')[[1]][9], '_')[[1]][[2]])))
    #file_list <- file_list[chunk <= 15]
    
    out <- do.call('c', lapply(file_list, readRDS))
    if (length(out) == 0) next()
    
    # Remove observations with errors
    check <- sapply(out, is.null)
    out <- out[!check]
    
    check <- sapply(out, function(x) inherits(x, 'error'))
    out <- out[!check]
    
    check <- sapply(out, anyNA, recursive=TRUE)
    out <- out[!check]
    
    check <- sapply(out, function(x) any(sapply(x, function(x) 
      ifelse(is.numeric(x), any(x > 1e5), FALSE))))
    out <- out[!check]
    
    check <- sapply(out, function(z) any(sapply(z, function(x) any(is.na(x)))))
    out <- out[!check]
    
    mc <- length(out)
    if (mc == 0) next()
    
    # Extract simualtion details
    n <- out[[1]]$n
    alpha <- out[[1]]$alpha
    design <- out[[1]]$design
    g1 <- out[[1]]$g1
    g2 <- out[[1]]$g2
    
    # Estimated parameters
    all_par <- do.call(rbind, lapply(out, '[[', 'par_fit'))
    all_par_qr <- do.call(rbind, lapply(out, '[[', 'par_fit_qr'))
    
    # Monte-Carlo covariance
    cov_mc <- cov(all_par) * n
    cov_mc_qr <- cov(all_par_qr) * n
    
    # Covariance estimators
    cov_est_names <- c('cov_iid_ind', 'cov_nid_scl_N', 'cov_nid_scl_sp', 'cov_boot')
    
    # Average Frobenius Norm
    avg_fnorm <- sapply(cov_est_names, function(x) 
      mean(sapply(lapply(out, "[[", x), function(z) fnorm(z * n - cov_mc))))
    
    # Average of the estimated covariances
    cov_est <- lapply(cov_est_names, function(x) 
      Reduce('+', lapply(out, '[[', x)) / length(out) * n)
    names(cov_est) <- cov_est_names
    
    cov_est_qr_names <- c('cov_qr_iid', 'cov_qr_nid')
    cov_est_qr <- lapply(cov_est_qr_names, 
                         function(x) Reduce('+', lapply(out, '[[', x)) / length(out) * n)
    names(cov_est_qr) <- cov_est_qr_names
    
    # Relative standard errors
    relse <- t(sapply(cov_est, function(cov) diag(cov_mc)^0.5 / diag(cov)^0.5))
    relse_qr <- t(sapply(cov_est_qr, function(cov) diag(cov_mc_qr)^0.5 / diag(cov)^0.5))
    
    # Times
    t_fit_mean <- mean(sapply(out, '[[', 't_fit'))
    t_fit_var <- var(sapply(out, '[[', 't_fit'))
    
    # True covariance
    true_cov_file <- paste0(sim_path, 'true_covariance/design_', design, '_g1_', g1, '_g2_', g2, '.rds')
    if (file.exists(true_cov_file) & (file.size(true_cov_file) > 0)) {
      df <- readRDS(true_cov_file)
      cov_true <- df$cov
      cov_true_qr <- df$cov_qr
    } else {
      cov_true <- matrix(NA, nrow(cov_est[[1]]), ncol(cov_est[[1]]))  
      cov_true_qr <- matrix(NA, nrow(cov_est_qr[[1]]), ncol(cov_est_qr[[1]]))  
    }
    
    # Stack everything
    results <- c(list(
      mc          = mc,
      n           = n,
      alpha       = alpha,
      design      = design,
      g1          = g1,
      g2          = g2,
      t_fit_mean  = t_fit_mean,
      t_fit_var   = t_fit_var,
      removed     = 1 - length(out) / mc,
      par_true    = out[[1]]$par_true,
      par_true_qr = out[[1]]$par_true[1:length(out[[1]]$par_fit_qr)],
      par_avg     = colMeans(all_par),
      par_avg_qr  = colMeans(all_par_qr),
      cov_true    = cov_true,
      cov_true_qr = cov_true_qr,
      cov_mc      = cov_mc,
      cov_mc_qr   = cov_mc_qr,
      relse       = relse,
      all_par     = all_par,
      avg_fnorm   = avg_fnorm
    ), cov_est, cov_est_qr)
    
    # Add to the big list
    all_results <- c(all_results, list(results))
  }
  
  all_results
}

get_identifier <- function(all_results) {
  id <- do.call('rbind', lapply(all_results, function(x)
    data.frame(mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
               stringsAsFactors = FALSE)))
  id
}

get_mse <- function(all_results) {
  mse <- do.call('rbind', lapply(all_results, function(x) 
    data.frame(mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
               mse=sum(t((x$par_true - x$par_avg)^2 + diag(x$cov_mc) / x$n)),
               stringsAsFactors = FALSE)))
  mse
}

get_mse_per_parameter <- function(all_results) {
  mse_per_parameter <- do.call('rbind', lapply(all_results, function(x) 
    data.frame(mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
               par=seq_len(nrow(x$cov_mc)),
               par_type=rep(c('Q', 'ES'), each=nrow(x$cov_mc) / 2),
               mse=(x$par_true - x$par_avg)^2 + diag(x$cov_mc) / x$n,
               stringsAsFactors = FALSE)))
  mse_per_parameter
}

get_frob_norm <- function(all_results) {
  frob_norm <- do.call('rbind', lapply(all_results, function(x) 
    data.frame(mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
               cov_estimator=names(x$avg_fnorm), norm=x$avg_fnorm,
               stringsAsFactors = FALSE)))
  
  frob_norm <- frob_norm %>% 
    group_by(design, g1, g2) %>%
    as.data.frame()
  
  frob_norm
}

get_relse <- function(all_results) {
  relse <- do.call('rbind', lapply(all_results, function(x) {
    data.frame(
      mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
      melt(x$relse, value.name = 'RELSE', varnames = c('cov_estimator', 'par'))
    )
  }))
  relse
}

get_estimation_time <- function(all_results) {
  estimation_time <- do.call('rbind', lapply(all_results, function(x) {
    data.frame(
      mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
      time=x$t_fit_mean
    )
  }))
  estimation_time
}

# Main part ---------------------------------------------------------------

# Folders with DGP's
sim_path <- '/mnt/Cluster/Repos/esreg_theory_paper/data/'
#sim_path <- '/mnt/Cluster/Repos_old/esreg_theory/data/'

all_results <- collect_all_results(sim_path)

# Get identifier
id <- get_identifier(all_results)

# Extract various tables
mse <- get_mse(all_results)
mse_per_parameter <- get_mse_per_parameter(all_results)
frob_norm <- get_frob_norm(all_results)
relse <- get_relse(all_results)
estimation_time <- get_estimation_time(all_results)

# Write tables to disk
write.csv(mse, file='../../Plots/Monte_Carlo/input/mse')
write.csv(mse_per_parameter, file='../../Plots/Monte_Carlo/input/mse_per_parameter')
write.csv(frob_norm, file='../../Plots/Monte_Carlo/input/norm')
write.csv(relse, file='../../Plots/Monte_Carlo/input/relse')
write.csv(estimation_time, file='../../Plots/Monte_Carlo/input/estimation_time')












# Comparison of MC and True Covariance ------------------------------------

cov_comp <- do.call('rbind', lapply(all_results, function(x) 
  data.frame(mc=x$mc, n=x$n, design=x$design, g1=x$g1, g2=x$g2,
             cov_true=fnorm(x$cov_true), cov_mc=fnorm(x$cov_mc),
             stringsAsFactors = FALSE)))

cov_comp %>% filter(n == 5000)


mse %>% 
  select(-mc) %>% 
  mutate(mse = sprintf("%0.3f", mse)) %>% 
  spread(g2, mse) %>% 
  arrange(design, n, g1)

mse_full %>% 
  select(-c(mc, par_type)) %>% 
  mutate(mse = sprintf("%0.3f", mse)) %>% 
  filter(design==1) %>% 
  spread(par, mse) %>% 
  arrange(g2, n)


mse_per_parameter %>% 
  select(-c(mc, par)) %>% 
  filter(design==1) %>% 
  group_by(n, design, g1, g2, par_type) %>% 
  summarise(mean_mse=mean(mse)) %>% 
  ungroup() %>% 
  mutate(mean_mse = sprintf("%0.3f", mean_mse)) %>% 
  spread(par_type, mean_mse) %>% 
  arrange(design, n, g1) %>% 
  as.data.frame()


estimation_time %>% 
  select(-c(mc)) %>% 
  mutate(time = sprintf("%0.3f", time)) %>% 
  spread(g2, time) 

frob_norm %>% 
  filter(design == 1, g2 == 1, n == 5000) %>% 
  spread(cov_estimator, norm) %>% 
  arrange(n)


mse_full %>% 
  filter(design==3, n==2000) %>% 
  select(-c(mc, par_type, design, n)) %>% 
  spread(g2, mse) 


# Asymptotic covariances ----------------------------------------------------------------------

replace_g1 <- function(g1) {
  switch(g1,
         '$G_1(z) = z$',
         '$G_1(z) = 0$',
         'Quantile Regression'
  )
}

replace_g2 <- function(g2) {
  switch(g2,
         '$\\mathcal{G}_2(z) = -\\log(-z)$',
         '$\\mathcal{G}_2(z) = -\\sqrt{-z}$',
         '$\\mathcal{G}_2(z) = -1/z$',
         '$\\mathcal{G}_2(z) = \\log(1 + \\exp(z))$',
         '$\\mathcal{G}_2(z) = \\exp(z)$',
         'Quantile Regression'
  )
}
for (g1 in c(1, 2)) {
  nn <- do.call('rbind', lapply(all_results, function(x) {
    if ((x$n == 250) & (x$g1 == g1)) {
      k <- nrow(x$cov_true) / 2
      xx <- data.frame(design=x$design, g1=x$g1, g2=x$g2, type=c(1:3),
                       norm=c(fnorm(x$cov_true[1:k, 1:k]),
                              fnorm(x$cov_true[(k+1):(2*k), (k+1):(2*k)]),
                              fnorm(x$cov_true)),
                       stringsAsFactors = FALSE)
      if (x$g2 == 1) {
        xx <- rbind(xx, data.frame(design=x$design, g1=g1, g2=6, type=1:3,
                                   norm=c(fnorm(x$cov_true_qr), -1, -1)))
      }
      
      xx
    }
  }))
  
  
  tab <- lapply(1:4, function(tt) t(nn %>% filter(design == tt) %>% spread(g2, norm) %>% select(-c(type, g1, design))))
  tab <- cbind(tab[[1]], NA, tab[[2]], NA, tab[[3]], NA, tab[[4]])
  
  rownames(tab) <- sapply(rownames(tab), function(x) replace_g2(as.numeric(x)))
  
  tab <- apply(tab, c(1, 2), function(x) {
    if (is.na(x)) {
      ""
    } else if (x == -1) {
      "--"
    } else if (x %% 1 == 0) {
      sprintf('%.0f', x)
    } else {
      sprintf('%.1f', x)
    }
  })
  
  xtab <- xtable(tab)
  print(xtab, #file = paste0('../../Tables/Monte_Carlo/true_variances_design_', dgp),
        sanitize.text.function = function(x) {x}, booktabs = TRUE,
        comment = FALSE, only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
}




# New Frob Norm -----------------------------------------------------------

