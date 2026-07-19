set.seed(42)
# Setup
n <- 200
k <- 1

beta_true <- c(1, 2)
sigma_eps <- 1.0
Omega <- matrix(0.25, nrow = 1, ncol = 1) 

B <- 1000
ols_estimates <- matrix(0, nrow = B, ncol = k+1)
mm_estimates <- matrix(0, nrow = B, ncol = k+1)

for (b in 1:B) {
  xi <- matrix(rnorm(n * k, mean = 0, sd = 1), nrow = n, ncol = k)
  vi <- matrix(rnorm(n * k, mean = 0, sd = sqrt(Omega)), nrow = n, ncol = k)
  xi_obs <- xi + vi
  eps <- rnorm(n, mean = 0, sd = sigma_eps)

  y <- cbind(1, xi) %*% beta_true + eps
  
  
  # OLS: (X'X)^-1 X'y
  xtx <- t(cbind(1, xi_obs)) %*% cbind(1, xi_obs)
  beta_ols <- solve(xtx) %*% t(cbind(1, xi_obs)) %*% y
  ols_estimates[b, ] <- beta_ols
  
  # MM: (X'X - n*Omega)^-1 X'y
  Omega_tilde <- matrix(0, nrow=k+1, ncol=k+1)
  Omega_tilde[2:(k+1), 2:(k+1)] <- Omega
  XtX_mm <- t(cbind(1, xi_obs)) %*% cbind(1, xi_obs) - n * Omega_tilde
  beta_mm <- solve(XtX_mm, t(cbind(1, xi_obs)) %*% y)
  mm_estimates[b, ] <- beta_mm
}


compute_metrics <- function(estimates, beta_true) {
  mean_est <- colMeans(estimates)
  bias <- mean_est - beta_true
  variance <- apply(estimates, 2, var) 
  mse <- bias^2 + variance
  list(mean=mean_est, bias=bias, var=variance, mse=mse)
}

ols_res <- compute_metrics(ols_estimates, beta_true)
mm_res  <- compute_metrics(mm_estimates, beta_true)


# result
param_names <- c("Intercept", paste0("X", 1:k))
results_df <- data.frame(
  Parameter = param_names,
  True_Beta = beta_true,
  # OLS 
  OLS_Mean = ols_res$mean,
  OLS_Bias = ols_res$bias,
  OLS_Var  = ols_res$var,
  OLS_MSE  = ols_res$mse,
  # MM 
  MM_Mean = mm_res$mean,
  MM_Bias = mm_res$bias,
  MM_Var  = mm_res$var,
  MM_MSE  = mm_res$mse
)
print(results_df)


out_dir <- "/home/u7080475/simulation/results2"
write.csv(results_df, file = file.path(out_dir, "res_olsmm_k=1.csv"), row.names = FALSE)
