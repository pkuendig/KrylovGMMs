################################################################################
# Runtime analysis for parameter estimation for a balanced random effects design
# (m_1=m_2) and Bernoulli likelihoods.
################################################################################

library(gpboost)

# options( krylovgmms.methods = c( "krylov",  "mixedmodels", "inla" ) )

.script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
.file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!nzchar(.script_path) && length(.file_args)) .script_path <- sub("^--file=", "", .file_args[1])
if (nzchar(.script_path)) setwd(dirname(normalizePath(.script_path)))
source("./../data/simulated/gen_data.R")
source("./../experiment_helpers.R")
if (method_enabled("lme4")) library(lme4)
if (method_enabled("glmmtmb")) library(glmmTMB)
set.seed(1)

M <- c(1000, 2000, 3000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)

sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2_1,sigma2_2)
init_cov_pars <- c(0.5,0.5)
num_covariates <- 5
init_betas <- c(0.0062, rep(0, num_covariates))

res_cols <- c("m","method", "time_estimation")
results <- data.frame(matrix(nrow=length(M)*6, ncol = length(res_cols)))
colnames(results) <- res_cols
log_time_vs_m_start(basename(.script_path), "Bernoulli", "balanced")

i <- 1
for(m in 1:length(M)){
  log_time_vs_m_iteration(M[m], "Bernoulli", "balanced")
  ###Generate data##############################################################
  set.seed(m)
  the_data <- make_data(n=M[m]*10, 
                        m1=M[m]/2,
                        sigma2=sigma2_1, #for signal-to-noise ratio
                        sigma2_1=sigma2_1,
                        sigma2_2=sigma2_2,
                        randef="Two_randomly_crossed_random_effects", 
                        likelihood="bernoulli_logit", 
                        has_F=TRUE,
                        num_covariates=num_covariates)

  ###Cholesky###################################################################
  results$m[i] <- M[m]
  results$method[i] <- "Cholesky (GPBoost)"

  if(method_enabled("cholesky") && M[m] <= 20000){
    log_time_vs_m_method(M[m], "cholesky")
    chol_model <- GPModel(group_data = the_data$group_data,
                          likelihood="bernoulli_logit",
                          matrix_inversion_method = "cholesky")

    chol_model$set_optim_params(params = list(maxit=1000,
                                              init_cov_pars=init_cov_pars,
                                              init_coef=init_betas))

    results$time_estimation[i] <- system.time(chol_model$fit(y=the_data$y, X=the_data$X))[3]
  }
  
  i <- i + 1
  
  ###Krylov#####################################################################
  results$m[i] <- M[m]
  results$method[i] <- "Krylov (GPBoost)"

  if (method_enabled("krylov")) {
  log_time_vs_m_method(M[m], "krylov")
  it_model <- GPModel(group_data = the_data$group_data,
                      likelihood="bernoulli_logit",
                      matrix_inversion_method = "iterative")
  
  it_model$set_optim_params(params = list(maxit=1000,
                                         init_cov_pars=init_cov_pars,
                                         init_coef=init_betas,
                                         seed_rand_vec_trace=1,
                                         cg_preconditioner_type="symmetric_successive_over_relaxation"))
  
  results$time_estimation[i] <- system.time(it_model$fit(y=the_data$y, X=the_data$X))[3]
  }
  
  i <- i + 1
  
  ###lme4#######################################################################
  results$m[i] <- M[m]
  results$method[i] <- "lme4"

  if(method_enabled("lme4") && M[m] <= 3000){
    log_time_vs_m_method(M[m], "lme4")
    formula = as.formula(paste0("y ~ -1 + ",paste0(colnames(the_data$X), collapse = ' + ')," + ",
                                paste0("(1|",colnames(the_data$group_data),")", collapse = ' + ')))
    results$time_estimation[i] <- system.time(lme4_model <- glmer(formula, family=binomial, data = data.frame(y = the_data$y, cbind(the_data$X, the_data$group_data)), 
                                                                  start=list(theta = lme4_start_theta(init_cov_pars, the_data$group_data, "bernoulli"),
                                                                             fixef = init_betas)))[3]
  }
    
  i <- i + 1
  
  ###glmmTMB####################################################################
  results$m[i] <- M[m]
  results$method[i] <- "glmmTMB"

  if(method_enabled("glmmtmb") && M[m] <= 10000){
    log_time_vs_m_method(M[m], "glmmtmb")
    ##Estimation
    formula = as.formula(paste0("y ~ -1 + ",paste0(colnames(the_data$X), collapse = ' + ')," + ",
                                paste0("(1|",colnames(the_data$group_data),")", collapse = ' + ')))
    results$time_estimation[i] <- system.time(glmmTMB_model <- glmmTMB(formula, family=binomial, data=data.frame(y = the_data$y, cbind(the_data$X, the_data$group_data)),
                                                                       start=list(theta = glmmtmb_start_theta(init_cov_pars, "bernoulli"),
                                                                                  beta = init_betas)))[3]
  }
    
  i <- i + 1

  ###MixedModels.jl#############################################################
  results$m[i] <- M[m]
  results$method[i] <- "MixedModels.jl"
  if (method_enabled("mixedmodels") && M[m] <= 20000) {
    log_time_vs_m_method(M[m], "mixedmodels")
    fit <- fit_additional_method("mixedmodels", the_data$y, the_data$X,
                                 the_data$group_data, "bernoulli",
                                 init_cov_pars, init_betas)
    results$time_estimation[i] <- fit$time_estimation
  }
  i <- i + 1

  ###R-INLA (empirical Bayes)###################################################
  results$m[i] <- M[m]
  results$method[i] <- "R-INLA (EB)"
  if (method_enabled("inla") && M[m] <= 5000) {
    log_time_vs_m_method(M[m], "inla")
    fit <- fit_additional_method("inla", the_data$y, the_data$X,
                                 the_data$group_data, "bernoulli",
                                 init_cov_pars, init_betas)
    results$time_estimation[i] <- fit$time_estimation
  }
  i <- i + 1
  log_time_vs_m_summary(M[m], results[(i - 6):(i - 1), ])
  saveRDS(results, "./../results/time_vs_m_balanced_Bernoulli.rds")
  gc()
}
