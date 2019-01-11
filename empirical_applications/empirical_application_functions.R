murphy <- function(r, q, e, alpha, v, return_mean=TRUE) {
  x <- (v <= e) * (1/alpha * (r <= q) * (q - r) - (q - v)) + (v <= r) * (r - v)
  if (return_mean) {
    mean(x)
  } else {
    x
  }
}

murpy_diff <- function(r, q1, e1, q2, e2, alpha, v) {
  m1 <- murphy(r=r, q=q1, e=e1, alpha=alpha, v=v, return_mean = FALSE)
  m2 <- murphy(r=r, q=q2, e=e2, alpha=alpha, v=v, return_mean = FALSE)
  m <- m1 - m2
  c(mean(m), sd(m) / sqrt(length(m)))
}

load_data <- function(file, symbol) {
  df <- read.csv(unz(file, "oxfordmanrealizedvolatilityindices.csv"), stringsAsFactors = FALSE)
  df_ <- df %>% dplyr::filter(Symbol == symbol) %>% dplyr::select(c('X', 'close_price', 'rv5'))
  y <- xts(diff.xts(log(df_[,2])), as.Date(df_[,1]))[-1] * 100
  x <- xts(sqrt(df_[,3]), as.Date(df_[,1]))[-1] * 100
  list(y = y, x = x)
}

