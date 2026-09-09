library(gpboost)

set.seed(1)

options( krylovgmms.methods = c("krylov") )

###Import data##################################################################
.source_paths <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else frame$ofile
}, character(1))
.script_path <- if (any(nzchar(.source_paths))) tail(.source_paths[nzchar(.source_paths)], 1) else ""
.file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!nzchar(.script_path) && length(.file_args)) .script_path <- sub("^--file=", "", .file_args[1])
if (!nzchar(.script_path)) .script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
if (nzchar(.script_path)) setwd(dirname(normalizePath(.script_path)))
source("./../experiment_helpers.R")
if (method_enabled("lme4")) library(lme4)
if (method_enabled("glmmtmb")) library(glmmTMB)
my_data <- read.csv("./../data/real_world/Amazon_employee_access.csv")

###Prepare data#################################################################
colnames(my_data) <- c("approved", "resource", "manager", "role_cat1", "role_cat2", 
                       "role_department", "role_title", "role_description", "role_family", 
                       "role_code")

#Categorical variables
cat_cols <- c("resource", "manager", "role_cat1", "role_cat2", 
              "role_department", "role_title", "role_description", "role_family", 
              "role_code")
group_data <- as.matrix(my_data[,cat_cols])

#Predictor variables
X <- as.matrix(rep(1, nrow(group_data)))
colnames(X)[1] <- "Intercept"

#Response
Y <- as.numeric(my_data$approved)

###Store data info##############################################################
data_cols <- c("ds_name", "n", "p", "K", "response_var", "class_imbalance", "cat_var", "nr_levels", "nnz_ZtZ")
data_info <- data.frame(matrix(nrow=length(cat_cols), ncol = length(data_cols)))
colnames(data_info) <- data_cols

data_info$ds_name[1] <- "amazon_employee_access"
data_info$n[1] <- length(Y)
data_info$p[1] <- ncol(X) - 1
data_info$K[1] <- ncol(group_data)
data_info$response_var[1] <- "approval/rejection"
data_info$class_imbalance[1] <- sum(Y)/length(Y)

for(i in 1:ncol(group_data)){
  data_info$cat_var[i] <- cat_cols[i]
  data_info$nr_levels[i] <- length(unique(group_data[,cat_cols[i]]))
}

# Infos on Z^TZ
data_info$nnz_ZtZ[1] <- count_group_ztz_nnz(group_data)
#image(ZtZ, xlab="", ylab="", sub ="", main="employee_access")

###Data frames##################################################################
res_cols <- c("method", "time_estimation", "num_optim_iter", "nll_optimum")
results <- data.frame(matrix(nrow=6, ncol = length(res_cols)))
colnames(results) <- res_cols
results$method <- c("Cholesky (GPBoost)", "Krylov (GPBoost)", "glmmTMB", "lme4",
                    "MixedModels.jl", "R-INLA (EB)")
log_real_world_start("amazon_employee_access")

cov_cols <- c(paste0("sigma2_", 1:length(cat_cols)))
cov_results <- data.frame(matrix(nrow=6, ncol = length(cov_cols)))
colnames(cov_results) <- cov_cols

beta_cols <- paste0("beta_", 0:(ncol(X)-1))
beta_results <- data.frame(matrix(nrow=6, ncol = length(beta_cols)))
colnames(beta_results) <- beta_cols

###Estimation###################################################################
i <- 1

# We use the same initial values for all libraries.
# The following are internal default initial values used by GPBoost
init_cov_pars <- rep(0.111111, ncol(group_data))
init_betas <- c(2.78958)

###Cholesky#####################################################################
if (method_enabled("cholesky")) {
  log_real_world_method("cholesky")
  chol_model <- GPModel(group_data = group_data,
                        likelihood="bernoulli_logit",
                        matrix_inversion_method = "cholesky")
  
  chol_model$set_optim_params(params = list(maxit=1000,
                                            trace=TRUE,
                                            init_cov_pars=init_cov_pars,
                                            init_coef=init_betas))
  
  results$method[i] <- "Cholesky (GPBoost)"
  results$time_estimation[i] <- system.time(chol_model$fit(y=Y, X=X))[3]
  beta_results[i,] <- chol_model$get_coef()
  cov_results[i,] <- chol_model$get_cov_pars()
  results$nll_optimum[i] <- chol_model$get_current_neg_log_likelihood()
  results$num_optim_iter[i] <- chol_model$get_num_optim_iter()
}
log_real_world_result("cholesky", results$time_estimation[i])

i <- i + 1

###Krylov#######################################################################
if (method_enabled("krylov")) {
  log_real_world_method("krylov")
  it_model <- GPModel(group_data = group_data,
                      likelihood="bernoulli_logit",
                      matrix_inversion_method = "iterative")
  
  it_model$set_optim_params(params = list(maxit=1000,
                                          trace=TRUE,
                                          init_cov_pars=init_cov_pars,
                                          init_coef=init_betas,
                                          seed_rand_vec_trace=50,
                                          cg_preconditioner_type="symmetric_successive_over_relaxation"))
  
  results$method[i] <- "Krylov (GPBoost)"
  results$time_estimation[i] <- system.time(it_model$fit(y=Y, X=X))[3]
  beta_results[i,] <- it_model$get_coef()
  cov_results[i,] <- it_model$get_cov_pars()
  results$nll_optimum[i] <- it_model$get_current_neg_log_likelihood()
  results$num_optim_iter[i] <- it_model$get_num_optim_iter()
}
log_real_world_result("krylov", results$time_estimation[i])

i <- i + 1

###glmmTMB####################################################################
if (method_enabled("glmmtmb")) {
  log_real_world_method("glmmtmb")
  formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_cols,")", collapse = ' + '))) 
  results$time_estimation[i] <- system.time(glmmTMB_model <- glmmTMB(formula, family=binomial, data = data.frame(y = Y, cbind(X, group_data)), 
                                                                     start=list(theta = glmmtmb_start_theta(init_cov_pars, "bernoulli"),
                                                                                beta = init_betas)))[3]
  
  results$method[i] <- "glmmTMB"
  beta_results[i,] <- fixef(glmmTMB_model)$cond
  cov_results[i,] <- as.numeric(VarCorr(glmmTMB_model)$cond)
  results$nll_optimum[i] <- -as.numeric(summary(glmmTMB_model)$logLik)
  results$num_optim_iter[i] <- glmmTMB_model$fit$iterations
}
log_real_world_result("glmmtmb", results$time_estimation[i])

i <- i + 1

###lme4#########################################################################
if (method_enabled("lme4")) {
  log_real_world_method("lme4")
  try({
    formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_cols,")", collapse = ' + '))) 
    results$time_estimation[i] <- system.time(lme4_model <- glmer(formula, family=binomial, data = data.frame(y = Y, cbind(X, group_data)), 
                                                                  start=list(theta = lme4_start_theta(init_cov_pars, group_data, "bernoulli"),
                                                                             fixef = init_betas)))[3]
    
    results$method[i] <- "lme4"
    beta_results[i,] <- summary(lme4_model)$coefficients[,1]
    vr <- as.data.frame(VarCorr(lme4_model))
    vr <- vr[match(cat_cols, vr$grp),]
    cov_results[i,] <- vr$vcov
    results$nll_optimum[i] <- -as.numeric(summary(lme4_model)$logLik)
    results$num_optim_iter[i] <- lme4_model@optinfo$feval
  })
}
log_real_world_result("lme4", results$time_estimation[i])

i <- i + 1
additional <- append_additional_results(results, cov_results, beta_results, i,
                                        Y, X, group_data, "bernoulli",
                                        init_cov_pars, init_betas)
results <- additional$results
cov_results <- additional$cov_results
beta_results <- additional$beta_results
log_real_world_summary(results)

################################################################################
all_results <- cbind(results, cov_results, beta_results)
save_real_world_results(all_results, data_info, "./../results/amazon_employee_access.rds")
