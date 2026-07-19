library(MASS) # mvrnorm
set.seed(42)
#setup
n <- 200
p_z <- 2  
p_x <- 2  

# true beta:(Intercept, Z1, Z2, X1, X2)
beta_true <- c(1, 0, 0, 0.25, -0.4) 
sigma_eps <- 1.0
k_total <- length(beta_true)
# 
Sigma_X <- matrix(c(1.0, -0.5,
                    -0.5, 1.0), nrow=2, byrow=TRUE)

Omega <- matrix(c(0.15, 0.1,
                  0.1, 0.75), nrow=2, byrow=TRUE)

B <- 1000
# check pd
stopifnot(all(eigen(Sigma_X)$values > 0))
stopifnot(all(eigen(Omega)$values > 0))

# MM 
correction_matrix <- matrix(0, nrow=k_total, ncol=k_total)
start_idx <- 1 + p_z + 1
end_idx   <- k_total
correction_matrix[start_idx:end_idx, start_idx:end_idx] <- Omega
print(correction_matrix)

#
ols_est <- matrix(0, nrow=B, ncol=k_total)
mm_est  <- matrix(0, nrow=B, ncol=k_total)

for (b in 1:B) {
  Z_binary <- rbinom(n, size = 1, prob = 0.5)
  Z_i <- matrix(rnorm(n * (p_z-1)), nrow=n, ncol=p_z-1)
  X_true <- mvrnorm(n, mu=rep(0, p_x), Sigma=Sigma_X)
  W_i <- mvrnorm(n, mu=rep(0, p_x), Sigma=Omega)
  X_obs <- X_true + W_i
  ones <- rep(1, n)
  D_true <- cbind(ones, Z_binary, Z_i, X_true)
  eps <- rnorm(n, mean=0, sd=sigma_eps)
  y <- D_true %*% beta_true + eps
  D_obs <- cbind(ones, Z_binary, Z_i, X_obs)
  
  # OLS (D'D)^(-1) D'y
  XtX_obs <- t(D_obs) %*% D_obs
  beta_ols <- solve(XtX_obs, t(D_obs) %*% y)
  ols_est[b, ] <- as.vector(beta_ols)
  
  # MM (D'D - n * Correction)^(-1) D'y
  XtX_corrected <- XtX_obs - n * correction_matrix
  beta_mm <- solve(XtX_corrected, t(D_obs) %*% y)
  mm_est[b, ] <- as.vector(beta_mm)
}

#
compute_metrics <- function(estimates, beta_true) {
  mean_est <- colMeans(estimates)
  bias <- mean_est - beta_true
  variance <- apply(estimates, 2, var) 
  mse <- bias^2 + variance
  list(mean=mean_est, bias=bias, var=variance, mse=mse)
}

ols_res <- compute_metrics(ols_est, beta_true)
mm_res  <- compute_metrics(mm_est, beta_true)

# result
param_names <- c("Intercept", paste0("Z", 1:p_z), paste0("X", 1:p_x))
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
write.csv(results_df, file = file.path(out_dir, "res_olsmm_withz_k=2.csv"), row.names = FALSE)



