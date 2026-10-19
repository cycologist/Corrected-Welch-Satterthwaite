script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0L) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  NA_character_
}
script_dir <- if (!is.na(script_path)) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}
source(file.path(script_dir, "corrected_welch_satterthwaite.R"))

set.seed(20260715)

K_values <- c(2L, 4L, 6L, 8L, 10L, 20L, 30L, 40L)
nu_values <- c(1L, 2L, 3L, 4L, 5L, 10L)
nsim <- 10000L

simulate_one <- function(K, nu, nsim) {
  u <- matrix(rchisq(nsim * K, df = nu), nrow = nsim, ncol = K)
  s2 <- u / nu
  df_est <- apply(s2, 1L, function(x) cws(x, nu = rep(nu, K)))
  target <- K * nu

  data.frame(
    K = K,
    nu = nu,
    target_df = target,
    mean_df = mean(df_est),
    bias = mean(df_est) - target,
    relative_bias = (mean(df_est) - target) / target,
    mc_se = sd(df_est) / sqrt(nsim),
    rmse = sqrt(mean((df_est - target)^2))
  )
}

grid <- expand.grid(K = K_values, nu = nu_values)
results <- do.call(
  rbind,
  Map(simulate_one, grid$K, grid$nu, MoreArgs = list(nsim = nsim))
)

print(results, row.names = FALSE, digits = 4)
