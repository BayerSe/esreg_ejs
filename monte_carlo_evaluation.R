rm(list = ls())

library(pacman)
p_load(naturalsort, dplyr, reshape2, tidyr, xtable, reticulate)
use_virtualenv('/home/sebastian/.virtualenvs/esreg_theory/') 
source_python('../esreg_python_functions/main.py')

tbl_folder <- '../../Revision EJS/Paper/v1/tables/'
img_folder <- '../../Revision EJS/Paper/v1/plots/'
sim_folder <- '/home/sebastian/Downloads/esreg_theory_data/'
max_chunk <- 1000


# Functions ---------------------------------------------------------------

concat_results <- function(all_results, var) {
  as_data_frame(do.call('rbind', lapply(all_results, function(x) get(var, x))))
}

fnorm <- function(mat, subset = NULL) {
  if (!is.null(subset)) {
    mat <- mat[subset, subset]  
  } 
  sqrt(sum(mat[lower.tri(mat, diag=TRUE)]^2))
}

collect_all_results <- function(sim_folder) {
  dir_list <- list.files(paste0(sim_folder, 'monte_carlo'), full.names = TRUE)
  dir_list <- naturalsort(dir_list)
  
  all_results <- list()
  for (idx in 1:length(dir_list)) {
    print(basename(dir_list[idx]))
    
    # Load results
    file_list <- list.files(dir_list[idx], full.names = TRUE)
    file_list <- file_list[file.size(file_list) > 0]
    
    chunks <- sapply(basename(file_list), function(x) as.numeric(strsplit(strsplit(x, '_')[[1]][2], '.', fixed = TRUE)[[1]][1]))
    file_list <- file_list[chunks <= max_chunk]
    n_files <- length(file_list)
    
    out <- do.call('c', lapply(file_list, readRDS))
    if (length(out) == 0) next()
    
    # Remove observations with errors
    check <- sapply(out, function(x) inherits(x, 'error'))
    out <- out[!check]
    
    check <- sapply(out, anyNA, recursive=TRUE)
    out <- out[!check]
    
    check <- apply(abs(do.call(rbind, lapply(out, '[[', 'par_fit'))) >= 100, 1, any)
    out <- out[!check]
    
    mc <- length(out)
    if (mc == 0) next()
    
    # Extract simulation settings
    n <- out[[1]]$n
    alpha <- out[[1]]$alpha
    design <- out[[1]]$design
    g1 <- out[[1]]$g1
    g2 <- out[[1]]$g2
    
    # Estimated parameters
    all_par <- do.call(rbind, lapply(out, '[[', 'par_fit'))
    all_par_qr <- do.call(rbind, lapply(out, '[[', 'par_fit_qr'))
    par_avg <- colMeans(all_par)
    
    # Monte-Carlo covariance
    cov_mc <- cov(all_par) * n
    cov_mc_qr <- cov(all_par_qr) * n
    #cov_mc <- robustbase::covMcd(all_par, nsamp = 'deterministic')$cov * n
    #cov_mc_qr <- robustbase::covMcd(all_par_qr, nsamp = 'deterministic')$cov * n
    
    # Mean squared error per parameter
    mse <- tibble(design=design, n=n, g1=g1, g2=g2,
                  par=seq_len(nrow(cov_mc)),
                  par_type=rep(c('Q', 'ES'), each=nrow(cov_mc) / 2),
                  value=(out[[1]]$par_true - par_avg)^2 + diag(cov_mc) / n)
    
    # Average Frobenius Norm
    cov_est_names <- c('cov_iid_ind', 'cov_nid_scl_N', 'cov_nid_scl_sp', 'cov_boot')
    avg_fnorm <- sapply(cov_est_names, function(x) 
      mean(sapply(lapply(out, "[[", x), function(z) fnorm(z * n - cov_mc))))
    avg_fnorm <- tibble(design=design, n=n, g1=g1, g2=g2, 
                        cov=cov_est_names, value=avg_fnorm)
    
    # True covariance
    true_cov_file <- paste0(sim_folder, 'true_covariance/design_', design, 
                            '_g1_', g1, '_g2_', g2, '.rds')
    df <- readRDS(true_cov_file)
    cov_true <- df$cov
    cov_true_qr <- df$cov_qr
    
    # Norm of true covariance
    k <- nrow(cov_true) / 2
    fnorm_true_cov <- tibble(design=design, n=n, g1=g1, g2=g2, 
                             type=c('Q', 'ES', 'Full'),
                             value=c(fnorm(cov_true[1:k, 1:k]),
                                     fnorm(cov_true[(k+1):(2*k), (k+1):(2*k)]),
                                     fnorm(cov_true)))
    fnorm_true_cov <- rbind(fnorm_true_cov, 
                            tibble(design=design, n=n, g1=g1, g2=6,
                                   type=c('Q', 'ES', 'Full'), 
                                   value=c(fnorm(cov_true_qr), NA, NA)))
    
    
    # # Average of the estimated covariances
    cov_est <- lapply(cov_est_names, function(x) 
      Reduce('+', lapply(out, '[[', x)) / length(out) * n)
    names(cov_est) <- cov_est_names
    
    # General information
    avg_time <- mean(sapply(out, '[[', 't_fit'))
    general_info <- tibble(design=design, n=n, g1=g1, g2=g2,
                           mc0=n_files * 100, mc=mc, 
                           rel_diff_mc = mc / (n_files * 100), 
                           time=avg_time)
    # Stack everything
    results <- list(
      general_info   = general_info,
      mse            = mse,
      avg_fnorm      = avg_fnorm,
      fnorm_true_cov = fnorm_true_cov
    )
    
    # Add to the big list
    all_results <- c(all_results, list(results))
  }
  
  all_results
}


# Main part ---------------------------------------------------------------

all_results <- collect_all_results(sim_folder)


# Extract results ---------------------------------------------------------

general_info <- concat_results(all_results, 'general_info')
mse_per_parameter <- concat_results(all_results, 'mse')
mse <- mse_per_parameter %>% 
  group_by(design, n, g1, g2) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()
frob_norm <- concat_results(all_results, 'avg_fnorm')
frob_norm_true_cov <- concat_results(all_results, 'fnorm_true_cov')


# Make plots --------------------------------------------------------------

design <- 3
n <- 2000
plot_mse_decomposed(mse_per_parameter, design=design, sample_size=n,
                    file=paste0(img_folder, 'mse_decomposed_design_', design, '_n_', n, '.pdf'))

plot_mse(mse, file=paste0(img_folder, 'mse.pdf'))

for (i in 1:2) {
  plot_norm(frob_norm, g1=i, file=paste0(img_folder, 'norm_g1_', i, '.pdf'))
}

plot_estimation_time(general_info, paste0(img_folder, 'estimation_time.pdf'))


# Asymptotic covariances ----------------------------------------------------------------------

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

for (design_ in list(c(1, 2), c(3, 4))) {
  for (g1_ in c(1, 2)) {
    tabs <- lapply(design_, function(d_) {
      frob_norm_true_cov %>% 
        subset(design == d_ & n == 250 & g1 == g1_) %>% 
        distinct(g2, type, .keep_all=TRUE) %>% 
        spread(type, value)
    })
    tab <- cbind(tabs[[1]][,c('Q', 'ES', 'Full')], -1, tabs[[2]][,c('Q', 'ES', 'Full')])
    rownames(tab) <- sapply(1:nrow(tab), function(x) replace_g2(as.numeric(x)))
    
    tab <- apply(tab, c(1, 2), function(x) {
      if (is.na(x)) {'--'} else if (x == -1) {''} else {sprintf('%.1f', x)}
    })
    
    file <- paste0(tbl_folder, 'true_covariance_design_', 
                   design_[1], design_[2], '_g1_', g1_, '.txt')
    xtab <- xtable(tab)
    print(xtab, file = file,
          sanitize.text.function = function(x) {x}, booktabs = TRUE,
          comment = FALSE, only.contents = TRUE, include.colnames = FALSE, hline.after = NULL)
    
  }
}

