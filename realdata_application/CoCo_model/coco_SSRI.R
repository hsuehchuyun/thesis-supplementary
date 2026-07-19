library(MASS)       # mvrnorm
library(hdme)       # correctlasso / mus / corrected lasso
library(glmnet)     # lasso
library(matrixcalc) # is.positive.definite
library(knitr)      

data <- read.csv("/home/u7080475/SSRI/FullW4_piHRSD_selected.csv")
y <- data$piHRSDW4
n <- nrow(data)
p_z <- 10  
p_x <- 50
k_total <- p_z + p_x
# ones <- rep(1, n)
Z <- as.matrix(data[, 4:13])
X_obs <- as.matrix(data[, 14:63])
D_obs <- cbind(Z, X_obs) 
# Omega_X <- as.matrix(read.csv("/home/u7080475/SSRI/CoCo_model/covariance/diff_SSRI_EUR(omega).csv", header = TRUE, row.names = 1))
Omega_X <- as.matrix(read.csv("/home/u7080475/SSRI/CoCo_model/covariance/diff_SSRI_EAS(omega).csv", header = TRUE, row.names = 1))

# standardize
sd_D <- apply(D_obs, 2, sd)
D_std <- scale(D_obs, center = TRUE, scale = TRUE)
y_std <- scale(y, center = TRUE, scale = TRUE)
S_inv <- diag(1 / sd_D)

#　Correction Matrix (p_z + p_x) x (p_z + p_x)
correction_matrix <- matrix(0, nrow=k_total, ncol=k_total)
start_idx <- p_z + 1
end_idx   <- k_total
correction_matrix[start_idx:end_idx, start_idx:end_idx] <- Omega_X
correction_matrix_std <- S_inv %*% correction_matrix %*% S_inv
rownames(correction_matrix_std) <- colnames(D_std)
colnames(correction_matrix_std) <- colnames(D_std)

# lasso
# choose lambda (cv, 10-fold, mse)
set.seed(1) 
cv_fit <- cv.glmnet(D_std, y_std, alpha = 1, family = "gaussian", intercept = FALSE ,nfolds = 5)
lambda_chosen <- cv_fit$lambda.min
fit_lasso <- glmnet(D_std, y_std, alpha = 1, family = "gaussian",intercept = FALSE,lambda = lambda_chosen)
beta_lasso <- as.numeric(coef(fit_lasso, s = lambda_chosen))
beta_lasso <- as.numeric(coef(fit_lasso, s = lambda_chosen))[-1]
print(beta_lasso)


# corrected lasso
set.seed(1)
fit_cv_corr <- cv_corrected_lasso(W = D_std, y = as.numeric(y_std), sigmaUU = correction_matrix_std, family = "gaussian",n_folds = 5)
print(fit_cv_corr)
# cv find best radius
print_text <- capture.output(print(fit_cv_corr))
target_line <- grep("Regularization parameter at minimum loss is", print_text, value = TRUE)
best_radius_corr <- as.numeric(regmatches(target_line, regexpr("[0-9]+\\.[0-9]+", target_line)))

# best_radius_corr <- 0.5444644 #EUR
# best_radius_corr <- 0.5305215 #EAS
fit_corr <- corrected_lasso(W = D_std, y = as.numeric(y_std), sigmaUU = correction_matrix_std, family = "gaussian", radii = best_radius_corr)
fit_corr$betaCorr

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

set.seed(1)

fit_cv_coco <- cv_corrected_lasso(W = D_std, y = as.numeric(y_std), sigmaUU = sigmaUU_coco_matrix, family = "gaussian",n_folds = 5)
print(fit_cv_coco)
# cv find best radius
print_text <- capture.output(print(fit_cv_coco))
target_line <- grep("Regularization parameter at minimum loss is", print_text, value = TRUE)
best_radius_coco <- as.numeric(regmatches(target_line, regexpr("[0-9]+\\.[0-9]+", target_line)))

# best_radius_coco <- 0.5514358 #EUR
# best_radius_coco <- 0.5374929 #EAS
fit_coco <- corrected_lasso(W = D_std, y = as.numeric(y_std), sigmaUU = sigmaUU_coco_matrix, family = "gaussian", radii = best_radius_coco)
fit_coco$betaCorr



# results
beta_lasso <- as.numeric(coef(fit_lasso, s = lambda_chosen))[-1]
beta_corr <- as.numeric(fit_corr$betaCorr)
beta_coco <- as.numeric(fit_coco$betaCorr)

var_names <- colnames(D_obs)
results_df <- data.frame(Variable = var_names, Lasso = beta_lasso, CorrectLasso = beta_corr, CoCoLasso = beta_coco)

output_path <- "/home/u7080475/SSRI/CoCo_model/result/coefficients_results(EAS).csv"
# output_path <- "/home/u7080475/SSRI/CoCo_model/result/coefficients_results(EUR).csv"
write.csv(results_df, file = output_path, row.names = FALSE)





