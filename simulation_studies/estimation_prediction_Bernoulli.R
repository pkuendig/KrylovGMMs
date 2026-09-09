################################################################################
# Repeated parameter estimation and prediction for Bernoulli likelihoods.
################################################################################

library(gpboost)
library(lme4)
library(glmmTMB)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")
source("./../experiment_helpers.R")
set.seed(1)

m1 <- 2000
n <- 2*(2*m1)*10
n_rep <- 100
sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2_1,sigma2_2)
init_cov_pars <- c(0.5,0.5)
num_covariates <- 5
init_betas <- c(0.0062, rep(0, num_covariates))

res_cols <- c("method",
              "RMSE_latent_mu", "mean_log_score",
              "sigma2_1", "sigma2_2",
              paste0("beta_", 0:num_covariates))
results <- data.frame(matrix(nrow=n_rep*3, ncol = length(res_cols)))
colnames(results) <- res_cols

i <- 1
for(r in 1:n_rep){
  ###Generate data##############################################################
  set.seed(r)
  the_data <- make_data(n=n, 
                        m1=m1,
                        sigma2=sigma2_1, #for signal-to-noise ratio
                        sigma2_1=sigma2_1,
                        sigma2_2=sigma2_2,
                        randef="Two_randomly_crossed_random_effects", 
                        likelihood="bernoulli_logit", 
                        has_F=TRUE,
                        num_covariates=num_covariates)
  
  #Sample train and test set
  i_train <- sample(n, n/2, replace = FALSE)
  i_test <- 1:n
  i_test <- i_test[!(i_test %in% i_train)]
  
  group_train <- the_data$group_data[i_train,]
  y_train <- the_data$y[i_train]
  eps_train <- the_data$eps[i_train]
  X_train <- the_data$X[i_train,]
  
  group_test <- the_data$group_data[i_test,]
  eps_test <- the_data$eps[i_test]
  X_test <- the_data$X[i_test,]
  
  ###Cholesky###################################################################
  
  ##Estimation
  chol_model <- GPModel(group_data = group_train,
                        likelihood="bernoulli_logit",
                        matrix_inversion_method = "cholesky")
  
  chol_model$set_optim_params(params = list(maxit=1000,
                                            init_cov_pars=init_cov_pars,
                                            init_coef=init_betas))
  
  chol_model$fit(y=y_train, X=X_train)

  results$method[i] <- "Cholesky (GPBoost)"
  results$sigma2_1[i] <- chol_model$get_cov_pars()[1]
  results$sigma2_2[i] <- chol_model$get_cov_pars()[2]
  betas_pred <- chol_model$get_coef()
  results[i,6:(num_covariates+6)] <- betas_pred
  
  ##Prediction
  pred_chol <- predict(chol_model,
                       group_data_pred=group_test,
                       X_pred = X_test,
                       predict_var = TRUE, 
                       predict_cov_mat = FALSE,
                       predict_response = FALSE)
  
  f_pred <- X_test %*% betas_pred
  eps_pred <- pred_chol$mu - f_pred
  
  results$mean_log_score[i] <- -mean(dnorm(eps_test, eps_pred, sqrt(pred_chol$var), log = TRUE))
  results$RMSE_latent_mu[i]  <- sqrt(mean((eps_pred-eps_test)^2))

  i <- i + 1
  
  ###Krylov#####################################################################

  ##Estimation
  it_model <- GPModel(group_data = group_train,
                      likelihood="bernoulli_logit",
                      matrix_inversion_method = "iterative")
  
  it_model$set_optim_params(params = list(maxit=1000,
                                         init_cov_pars=init_cov_pars,
                                         init_coef=init_betas,
                                         seed_rand_vec_trace=r,
                                         cg_preconditioner_type="symmetric_successive_over_relaxation"))
  
  it_model$fit(y=y_train, X=X_train)
  
  results$method[i] <- "Krylov (GPBoost)"
  results$sigma2_1[i] <-  it_model$get_cov_pars()[1]
  results$sigma2_2[i] <-  it_model$get_cov_pars()[2]
  betas_pred <- it_model$get_coef()
  results[i,6:(num_covariates+6)] <- betas_pred
  
  ##Prediction
  pred_it <- predict(it_model,
                     group_data_pred=group_test,
                     X_pred = X_test,
                     predict_var = TRUE, 
                     predict_cov_mat = FALSE,
                     predict_response = FALSE)
  
  f_pred <- X_test %*% betas_pred
  eps_pred <- pred_it$mu - f_pred
  
  results$mean_log_score[i] <- -mean(dnorm(eps_test, eps_pred, sqrt(pred_it$var), log = TRUE))
  results$RMSE_latent_mu[i]  <- sqrt(mean((eps_pred-eps_test)^2))

  i <- i + 1

  ###glmmTMB####################################################################
  
  ##Estimation  
  formula = as.formula(paste0("y ~ -1 + ",paste0(colnames(X_train), collapse = ' + ')," + ",
                              paste0("(1|",colnames(group_train),")", collapse = ' + ')))
  glmmTMB_model <- glmmTMB(formula, family=binomial, data=data.frame(y = y_train, cbind(X_train, group_train)),
                           start=list(theta = glmmtmb_start_theta(init_cov_pars, "bernoulli"),
                                      beta = init_betas))
  
  results$method[i] <- "glmmTMB"
  results$sigma2_1[i] <- as.numeric(VarCorr(glmmTMB_model)$cond)[1]
  results$sigma2_2[i] <- as.numeric(VarCorr(glmmTMB_model)$cond)[2]
  betas_pred <- fixef(glmmTMB_model)$cond
  results[i,6:(num_covariates+6)] <- betas_pred
  
  ##Prediction
  pred_glmmTMB <- predict(glmmTMB_model, type = "link", se.fit = FALSE, newdata = as.data.frame(cbind(X_test, group_test)))
  
  f_pred <- X_test %*% betas_pred
  eps_pred <- pred_glmmTMB - f_pred
  
  results$RMSE_latent_mu[i]  <- sqrt(mean((eps_pred-eps_test)^2)) 

  i <- i + 1
  
  saveRDS(results, "./estimation_prediction_Bernoulli.rds")
  gc()
}
