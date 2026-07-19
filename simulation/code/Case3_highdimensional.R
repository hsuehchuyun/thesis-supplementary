library(MASS)       # mvrnorm
library(hdme)       # correctlasso / mus / corrected lasso
library(glmnet)     # lasso
library(matrixcalc) # is.positive.definite
library(knitr)      

set.seed(42)
# 
n   <- 200
p_z <- 10  
p_x <- 50  
B   <- 1000  
sigma_eps <- 1.0

# true beta (Intercept) + 10 (Z) + 50 (X)
k_total <- 1 + p_z + p_x
beta_true <- rep(0, k_total)
beta_true[1] <- 1  #(Intercept)

# Z  Z1-Z3!=0, Z4-Z10 = 0
# beta_true[2:4] <- c(2,-1,0.5)
beta_true[2:4] <- rnorm(3) * 10 

# X  X1-X5!=0, X6-X50 = 0 
beta_true[12:16] <- rnorm(5) * 10 

# Sigma_X 
df_X <- p_x + 5
df_X <- p_x 
W_X  <- rWishart(1, df=df_X, Sigma=diag(p_x))[,,1]
D_inv <- diag(1 / sqrt(diag(W_X)))
Sigma_X <- D_inv %*% W_X %*% D_inv
Sigma_X <- (Sigma_X + t(Sigma_X)) / 2 
stopifnot(is.positive.definite(Sigma_X))

# Omega
var_val <- 0.5  
cov_val <- 0.1 
Omega_X <- matrix(cov_val, nrow = p_x, ncol = p_x)
diag(Omega_X) <- var_val
stopifnot(is.positive.definite(Omega_X))
# print(Omega_X)

# MM (Correction Matrix) (1 + p_z + p_x) x (1 + p_z + p_x)
correction_matrix <- matrix(0, nrow=k_total, ncol=k_total)

#
start_idx <- 1 + p_z + 1 
end_idx   <- k_total
correction_matrix[start_idx:end_idx, start_idx:end_idx] <- Omega_X

# save data
ols_est <- matrix(0, nrow=B, ncol=k_total)
mm_est  <- matrix(0, nrow=B, ncol=k_total)
lasso_est  <- matrix(0, nrow=B, ncol=k_total)
correctlasso_est  <- matrix(0, nrow=B, ncol=k_total)
cocolasso_est <- matrix(0, nrow=B, ncol=k_total)


for (b in 1:B) {
  set.seed(100 + b)
  Z_binary <- rbinom(n, size = 1, prob = 0.5)
  Z_i <- matrix(rnorm(n * (p_z-1)), nrow=n, ncol=p_z-1)
  X_true <- mvrnorm(n, mu=rep(0, p_x), Sigma=Sigma_X)
  V_i <- mvrnorm(n, mu=rep(0, p_x), Sigma=Omega_X)
  X_obs <- X_true + V_i
  
  # Y
  ones <- rep(1, n)
  D_true <- cbind(ones, Z_binary, Z_i, X_true) 
  eps <- rnorm(n, mean=0, sd=sigma_eps)
  y <- D_true %*% beta_true + eps
  
  #
  D_obs <- cbind(ones, Z_binary, Z_i, X_obs)   
  XtX_obs <- t(D_obs) %*% D_obs      # X'X
  Xty_obs <- t(D_obs) %*% y          # X'y
  # ---method---
  # OLS: (X'X)^-1 X'y
  beta_ols <- solve(XtX_obs, Xty_obs)
  ols_est[b, ] <- as.vector(beta_ols)
  
  # MM: (X'X - n*Omega)^-1 X'y
  XtX_corrected <- XtX_obs - n * correction_matrix
  beta_mm <- solve(XtX_corrected, Xty_obs)
  mm_est[b, ] <- as.vector(beta_mm)
  
  # LASSO
  # choose lamda
  # lambda_chosen <- 2   #
  Predictors <- cbind(Z_binary, Z_i, X_obs)
  cv_fit <- cv.glmnet(Predictors, y, alpha = 1, family = "gaussian" ,nfolds = 10) #cv
  lambda_chosen <- cv_fit$lambda.min #cv
  fit_lasso <- glmnet(Predictors, y, alpha = 1, family = "gaussian",lambda = lambda_chosen)
  beta_hat <- as.numeric(coef(fit_lasso, s = lambda_chosen))
  lasso_est[b, ] <- beta_hat
  
  # correctlasso
  # fit_corr <- corrected_lasso(D_obs, y, sigmaUU = correction_matrix, family = "gaussian")
  # # choose no. 6 radius
  # my_radius <- fit_corr$radii[6]

  fit_cv_corr <- cv_corrected_lasso(D_obs, y, sigmaUU = correction_matrix, family = "gaussian",n_folds = 5)
  print(fit_cv_corr)
  # cv find best radius
  print_text <- capture.output(print(fit_cv_corr))
  target_line <- grep("Regularization parameter at minimum loss is", print_text, value = TRUE)
  my_radius <- as.numeric(regmatches(target_line, regexpr("[0-9]+\\.[0-9]+", target_line)))

  fit_corr_single <- corrected_lasso(D_obs, y, sigmaUU = correction_matrix, family = "gaussian", radii = my_radius)
  correctlasso_est[b, ] <-  fit_corr_single$betaCorr
  
  # CoCoLasso
  Sigma_hat <- (t(X_obs) %*% X_obs) / n - Omega_X
  
  #ADMM
  admm_appendix_a <- function(S, mu = 3, epsilon = 1e-6, tol = 1e-4, maxiter = 1000) {
    #  S 對稱
    S <- (S + t(S)) / 2
    # 初始化
    B <- S # B0
    Lambda <- matrix(0, nrow(S), ncol(S)) # Lambda0
    A <- matrix(0, nrow(S), ncol(S)) # Ai+1 初始化
    
    # 迭代
    for (i in 1:maxiter) {
      B_old <- B
      # ----------------------------------------------------
      #step A (Ai+1=(Bi+Σ^+μΛi)ε)
      target_A <- B + S + mu * Lambda
      eig_A <- eigen(target_A, symmetric = TRUE) #特徵分解（Spectral decomposition）
      values_A_proj <- pmax(eig_A$values, epsilon) # 特徵值太小的，提升到至少是𝜀(修正成半正定)
      A <- eig_A$vectors %*% diag(values_A_proj) %*% t(eig_A$vectors) #重建矩陣
      
      # ----------------------------------------------------
      # Step B : Bi+1 = ... l_1(...)
      target_B <- A - S - mu * Lambda
      # 將矩陣向量化 (取下三角)
      vec_L_target_B <- target_B[lower.tri(target_B, diag = TRUE)] #vecl()
      # l_1(x, mu) -- L1 Ball 投影 (hdme套件裡的函示)
      l1_projection <- hdme:::project_onto_l1_ball(vec_L_target_B, mu) #l_1()
      # vecl(...) - l_1(...)
      vec_L_B_new <- vec_L_target_B - l1_projection
      # 重建矩陣 matl(...)
      B <- matrix(0, nrow(S), ncol(S))
      B[lower.tri(B, diag = TRUE)] <- vec_L_B_new
      B <- B + t(B)
      diag(B) <- diag(B) / 2 # 對角線元素在 t(B) 裡被加了兩次
      
      # ----------------------------------------------------
      # Step Lambda
      Lambda <- Lambda - (A - B - S) / mu
      
      # 檢查收斂
      if (norm(B - B_old, "F") < tol) {
        cat("Converged in", i, "iterations.\n")
        break
      }
    }
    return(A)
  }
  admmSigma <- admm_appendix_a(Sigma_hat)
  eigen(admmSigma)$values
  
  # ---- Step 3: 更新 sigmaUU ----
  sigmaUU_coco <- (t(X_obs) %*% X_obs) / n - admmSigma
  sigmaUU_coco_matrix <- matrix(0, nrow=k_total, ncol=k_total)
  sigmaUU_coco_matrix[start_idx:end_idx, start_idx:end_idx] <- sigmaUU_coco
  
  # ---- Step 4: 用 corrected_lasso() ----
  # fit_coco <- corrected_lasso(D_obs, y, sigmaUU = sigmaUU_coco_matrix, family = "gaussian")
  # my_radius <- fit_coco$radii[6]
  fit_cv_coco <- cv_corrected_lasso(D_obs, y, sigmaUU = sigmaUU_coco_matrix, family = "gaussian",n_folds = 5)
  print(fit_cv_coco)
  # cv find best radius
  print_text <- capture.output(print(fit_cv_coco))
  target_line <- grep("Regularization parameter at minimum loss is", print_text, value = TRUE)
  my_radius <- as.numeric(regmatches(target_line, regexpr("[0-9]+\\.[0-9]+", target_line)))
  fit_coco_single <- corrected_lasso(D_obs, y, sigmaUU = sigmaUU_coco_matrix,family = "gaussian", radii = my_radius)
  cocolasso_est[b, ] <- fit_coco_single$betaCorr
}

# 
compute_metrics <- function(estimates, beta_true) {
  mean_est <- colMeans(estimates)
  bias <- mean_est - beta_true
  var <- apply(estimates, 2, var)
  mse <- bias^2 + var
  l2norm <- sum(bias^2) 
  list(mean=mean_est, bias=bias, var=var, mse=mse, l2norm=l2norm)
}

ols_res_withz <- compute_metrics(ols_est, beta_true)
mm_res_withz  <- compute_metrics(mm_est, beta_true)
lasso_res_withz  <- compute_metrics(lasso_est, beta_true)
correctlasso_res_withz  <- compute_metrics(correctlasso_est, beta_true)
cocolasso_res_withz <- compute_metrics(cocolasso_est, beta_true)


# result
param_names <- c("Intercept", paste0("Z", 1:p_z), paste0("X", 1:p_x))

res_highdim_withz <- data.frame(
  Parameter = param_names,
  True_Beta = beta_true,
  OLS_Mean = ols_res_withz$mean,
  OLS_Bias = ols_res_withz$bias,
  OLS_Var = ols_res_withz$var,
  OLS_MSE  = ols_res_withz$mse,
  OLS_L2 = ols_res_withz$l2norm,
  MM_Mean = mm_res_withz$mean,
  MM_Bias = mm_res_withz$bias,
  MM_Var = mm_res_withz$var,
  MM_MSE  = mm_res_withz$mse,
  MM_L2 = mm_res_withz$l2norm,
  LASSO_Mean = lasso_res_withz$mean,
  LASSO_Bias = lasso_res_withz$bias,
  LASSO_Var = lasso_res_withz$var,
  LASSO_MSE  = lasso_res_withz$mse,
  LASSO_L2 = lasso_res_withz$l2norm,
  CorrectLASSO_Mean = correctlasso_res_withz$mean,
  CorrectLASSO_Bias = correctlasso_res_withz$bias,
  CorrectLASSO_Var = correctlasso_res_withz$var,
  CorrectLASSO_MSE  = correctlasso_res_withz$mse,
  CorrectLASSO_L2 = correctlasso_res_withz$l2norm,
  CoCoLASSO_Mean = cocolasso_res_withz$mean,
  CoCoLASSO_Bias = cocolasso_res_withz$bias,
  CoCoLASSO_Var = cocolasso_res_withz$var,
  CoCoLASSO_MSE  = cocolasso_res_withz$mse,
  CoCoLASSO_L2 = cocolasso_res_withz$l2norm
)
print(res_highdim_withz)

out_dir <- "/home/u7080475/simulation/results2"
write.csv(res_highdim_withz, file = file.path(out_dir, "res_highdim_withz.csv"), row.names = FALSE)

get_selection_metrics <- function(estimates, beta_true, tolerance = 1e-5) {
  est_coefs <- estimates[, 2:ncol(estimates)]
  true_coefs <- beta_true[2:length(beta_true)]
  B <- nrow(est_coefs)
  idx_true_nonzero <- which(abs(true_coefs) > tolerance) # 真正有作用的變數索引
  idx_true_zero    <- which(abs(true_coefs) <= tolerance) # 真正為 0 的變數索引
  
  total_non_zeros <- 0
  true_non_zeros  <- 0
  false_non_zeros <- 0
  
  for (b in 1:B) {
    beta_hat <- est_coefs[b, ]
    idx_est_nonzero <- which(abs(beta_hat) > tolerance)
    total_non_zeros <- total_non_zeros + length(idx_est_nonzero)
    true_non_zeros <- true_non_zeros + length(intersect(idx_est_nonzero, idx_true_nonzero))
    false_non_zeros <- false_non_zeros + length(intersect(idx_est_nonzero, idx_true_zero))
  }

  list(
    avg_nz = total_non_zeros / B,
    avg_true_nz = true_non_zeros / B,
    avg_false_nz = false_non_zeros / B
  )
}

sel_ols          <- get_selection_metrics(ols_est, beta_true)
sel_mm           <- get_selection_metrics(mm_est, beta_true)
sel_lasso        <- get_selection_metrics(lasso_est, beta_true)
sel_correctlasso <- get_selection_metrics(correctlasso_est, beta_true)
sel_cocolasso    <- get_selection_metrics(cocolasso_est, beta_true)

#
row_names <- c(
  "L2 norm",
  "Average number of non-zeros estimated",
  "Average number of true non-zeros estimated",
  "Average number of false non-zeros estimated"
)

# 
col_ols <- c(ols_res_withz$l2norm, sel_ols$avg_nz, sel_ols$avg_true_nz, sel_ols$avg_false_nz)
col_mm <- c(mm_res_withz$l2norm, sel_mm$avg_nz, sel_mm$avg_true_nz, sel_mm$avg_false_nz)
col_lasso <- c(lasso_res_withz$l2norm, sel_lasso$avg_nz, sel_lasso$avg_true_nz, sel_lasso$avg_false_nz)
col_correctlasso <- c(correctlasso_res_withz$l2norm, sel_correctlasso$avg_nz, sel_correctlasso$avg_true_nz, sel_correctlasso$avg_false_nz)
col_cocolasso <- c(cocolasso_res_withz$l2norm, sel_cocolasso$avg_nz, sel_cocolasso$avg_true_nz, sel_cocolasso$avg_false_nz)
# Data Frame
final_table <- data.frame(
  Metric = row_names,
  OLS = round(col_ols, 4),
  MM = round(col_mm, 4),
  Lasso = round(col_lasso, 4),
  Corrected_Lasso = round(col_correctlasso, 4),
  CoCoLasso = round(col_cocolasso, 4)
)

print(final_table)

out_dir <- "/home/u7080475/simulation/results2"
write.csv(final_table, file = file.path(out_dir, "variable_selection_withz.csv"), row.names = FALSE)


