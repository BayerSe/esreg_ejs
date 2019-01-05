rm(list = ls())

library(pacman)
p_load(esreg, foreach, doParallel, quantreg)

source("dgp_functions.R")

# Start the pool
nodename <- Sys.info()["nodename"]
if (grepl("uc1", nodename)) {
  cores <- as.numeric(Sys.getenv("SLURM_NPROCS"))
  registerDoParallel(cores=cores)
} else {
  registerDoSEQ() 
}

# Make path
path <- "data/monte_carlo/"
dir.create(path, showWarnings = FALSE, recursive = TRUE)

# Simulation setup
setup <- expand.grid(
  alpha   = 0.025,
  n       = c(250, 500, 1000, 2000, 5000),
  design  = c(1, 2, 3, 4),
  g1      = c(1, 2),
  g2      = c(1, 2, 3, 4, 5)
)

# Randomize the order of execution
setup <- setup[sample(nrow(setup)),]

# Number of Monte-Carlo repetitions and the chunks
mc <- 1:25000
save_frequency <- 100
all_chunks <- split(mc, ceiling(seq_along(mc) / save_frequency))

# Run the simulation: loop over all chunks and settings
for (chunk in 1:length(all_chunks)) {
  for (idx in 1:nrow(setup)) {
    st <- setup[idx,]
    
    # Set filename and check whether it already exists
    dir <- paste0(path, "design_", st$design, 
                  "_n_", st$n, "_g1_", st$g1, "_g2_", st$g2, "/")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    file <- paste0(dir, "chunk_", chunk, ".rds")
    if (file.exists(file)) next
    print(paste0("Computing ", file))
    tmp <- file.create(file)
    
    # Compute the current set of seeds
    out <- foreach(seed = all_chunks[[chunk]], .errorhandling = "pass") %dopar% {
      # Simulate from the dgp
      df <- dgp_fun(design = st$design, n = st$n, alpha = st$alpha, seed = seed)
      
      # Fit model and estimate the covariances
      fit <- esreg(df$y ~ df$xq - 1 | df$xe - 1, alpha = st$alpha, g1 = st$g1, g2 = st$g2)
      all_cov <- list(
        cov_iid_ind     = vcovA(fit, sparsity = "iid", cond_var = "ind"),
        cov_nid_scl_N   = vcovA(fit, sparsity = "nid", cond_var = "scl_N"),
        cov_nid_scl_sp  = vcovA(fit, sparsity = "nid", cond_var = "scl_sp"),
        cov_boot        = vcovB(fit, B = 2000, bootstrap_method = "iid")
      )
      
      # Quantile regression for comparison
      fit_qr <- rq(df$y ~ df$xq - 1, tau = st$alpha)
      all_cov_qr <- list(
        cov_qr_iid = summary.rq(fit_qr, covariance = TRUE, hs = TRUE, se = "iid")$cov,
        cov_qr_nid = summary.rq(fit_qr, covariance = TRUE, hs = TRUE, se = "nid")$cov
      )
      
      # Return
      c(n = st$n, alpha = st$alpha, design = st$design, 
        g1 = st$g1, g2 = st$g2, seed = seed, 
        t_fit      = fit$time, 
        par_true   = list(c(df$par_q, df$par_e)), 
        par_fit    =  list(unname(fit$coefficients)), 
        par_fit_qr =  list(unname(fit_qr$coefficients)),
        lapply(all_cov, unname), 
        lapply(all_cov_qr, unname))
    }
    
    # Store the results
    saveRDS(out, file)
  }
}
