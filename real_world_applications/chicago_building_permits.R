library(gpboost)

set.seed(1)

# "mixedmodels" crashes
options( krylovgmms.methods = c("cholesky", "krylov", "lme4", "glmmtmb", "inla") )

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
my_data <- read.csv(unz("./../data/real_world/chicago_building_permits.zip","chicago_building_permits.csv"))

my_data <- my_data[!is.na(my_data$REPORTED_COST),]
my_data <- my_data[my_data$REPORTED_COST<=1e7,]
my_data <- my_data[my_data$REPORTED_COST>1,]

###Prepare data#################################################################

#Categorical variables
colnames(my_data)[27] <- "contact_name"
colnames(my_data)[117] <- "latitude"
colnames(my_data)[118] <- "longitude"

cat_vars <- c("contact_name", "latitude", "longitude")
group_data_df <- my_data[,cat_vars]
group_data <- matrix(nrow = nrow(group_data_df), ncol=ncol(group_data_df))
colnames(group_data) <- colnames(group_data_df)
for(j in 1:dim(group_data)[2]) {
  group_data[,j] <- as.numeric(addNA(group_data_df[,j]))
}

#Predictor variables
pred_vars <- c("PROCESSING_TIME", "APPLICATION_START_DATE", "STREET.DIRECTION", "TOTAL_FEE")

my_data$PROCESSING_TIME[my_data$PROCESSING_TIME<0] <- 0
my_data$PROCESSING_TIME[is.na(my_data$PROCESSING_TIME)] <- mean(my_data$PROCESSING_TIME, na.rm=TRUE)

my_data$APPLICATION_START_DATE <- as.numeric(as.Date(my_data$APPLICATION_START_DATE))
my_data$APPLICATION_START_DATE[is.na(my_data$APPLICATION_START_DATE)] <- mean(my_data$APPLICATION_START_DATE, na.rm=TRUE)

my_data$STREET.DIRECTION[is.na(my_data$STREET.DIRECTION)] <- "W"
my_data$STREET.DIRECTION <- as.factor(my_data$STREET.DIRECTION)

my_data$TOTAL_FEE[my_data$TOTAL_FEE<1] <- 1
my_data$TOTAL_FEE <- log(my_data$TOTAL_FEE)

X <- model.matrix(~., my_data[,pred_vars])
colnames(X)[1] <- "Intercept"

#Response
Y <- log(my_data$REPORTED_COST)
rm(my_data)

###Store data info##############################################################
data_cols <- c("ds_name", "n", "p", "K", "response_var", "class_imbalance", "cat_var", "nr_levels", "nnz_ZtZ")
data_info <- data.frame(matrix(nrow=length(cat_vars), ncol = length(data_cols)))
colnames(data_info) <- data_cols

data_info$ds_name[1] <- "chicago_building_permits"
data_info$n[1] <- length(Y)
data_info$p[1] <- ncol(X) - 1
data_info$K[1] <- ncol(group_data)
data_info$response_var[1] <- "log(cost)"

for(i in 1:ncol(group_data)){
  data_info$cat_var[i] <- cat_vars[i]
  data_info$nr_levels[i] <- length(unique(group_data[,cat_vars[i]]))
}

# Infos on Z^TZ
data_info$nnz_ZtZ[1] <- count_group_ztz_nnz(group_data)
#image(ZtZ, xlab="", ylab="", sub ="", main="building_permits")

###Data frames##################################################################
res_cols <- c("method", "time_estimation", "num_optim_iter", "nll_optimum")
results <- data.frame(matrix(nrow=6, ncol = length(res_cols)))
colnames(results) <- res_cols
results$method <- c("Cholesky (GPBoost)", "Krylov (GPBoost)", "glmmTMB", "lme4",
                    "MixedModels.jl", "R-INLA (EB)")
log_real_world_start("chicago_building_permits")

cov_cols <- c("sigma", paste0("sigma2_", 1:length(cat_vars)))
cov_results <- data.frame(matrix(nrow=6, ncol = length(cov_cols)))
colnames(cov_results) <- cov_cols

beta_cols <- paste0("beta_", 0:(ncol(X)-1))
beta_results <- data.frame(matrix(nrow=6, ncol = length(beta_cols)))
colnames(beta_results) <- beta_cols

###Estimation###################################################################
i <- 1

# We use the same initial values for all libraries.
# The following are internal default initial values used by GPBoost
init_cov_pars <- c(2.81275, rep(0.937582, ncol(group_data)))
init_betas <- c(8.66558, rep(0, ncol(X)-1))

###Cholesky#####################################################################
if (method_enabled("cholesky")) {
  log_real_world_method("cholesky")
  chol_model <- GPModel(group_data = group_data,
                        likelihood="gaussian",
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
                      likelihood="gaussian",
                      matrix_inversion_method = "iterative")
  
  it_model$set_optim_params(params = list(maxit=1000,
                                          trace=TRUE,
                                          init_cov_pars=init_cov_pars,
                                          init_coef=init_betas,
                                          cg_preconditioner_type="zero_infill_incomplete_cholesky"))
  
  results$method[i] <- "Krylov (GPBoost)"
  results$time_estimation[i] <- system.time(it_model$fit(y=Y, X=X))[3]
  beta_results[i,] <- it_model$get_coef()
  cov_results[i,] <- it_model$get_cov_pars()
  results$nll_optimum[i] <- it_model$get_current_neg_log_likelihood()
  results$num_optim_iter[i] <- it_model$get_num_optim_iter()
}
log_real_world_result("krylov", results$time_estimation[i])

i <- i + 1

###glmmTMB######################################################################
if (method_enabled("glmmtmb")) {
  log_real_world_method("glmmtmb")
  formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
  results$time_estimation[i] <- system.time(glmmTMB_model <- glmmTMB(formula, family=gaussian(), data = data.frame(y = Y, cbind(X, group_data)), 
                                                                     start=list(theta = glmmtmb_start_theta(init_cov_pars, "gaussian"),
                                                                                betadisp = glmmtmb_start_betadisp(init_cov_pars),
                                                                                beta = init_betas)))[3]
  
  results$method[i] <- "glmmTMB"
  beta_results[i,] <- fixef(glmmTMB_model)$cond
  cov_results[i,] <- c((summary(glmmTMB_model)$sigma)^2, as.numeric(VarCorr(glmmTMB_model)$cond))
  results$nll_optimum[i] <- -as.numeric(summary(glmmTMB_model)$logLik)
  results$num_optim_iter[i] <- glmmTMB_model$fit$iterations
}
log_real_world_result("glmmtmb", results$time_estimation[i])

i <- i + 1

###lme4#########################################################################
if (method_enabled("lme4")) {
  log_real_world_method("lme4")
  try({
    formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
    results$time_estimation[i] <- system.time(lme4_model <- lmer(formula, REML = FALSE, data = data.frame(y = Y, cbind(X, group_data)), 
                                                                 start=list(theta = lme4_start_theta(init_cov_pars, group_data, "gaussian"))))[3]
    
    results$method[i] <- "lme4"
    beta_results[i,] <- summary(lme4_model)$coefficients[,1]
    vr <- as.data.frame(VarCorr(lme4_model))
    vr <- vr[match(c("Residual", cat_vars), vr$grp),]
    cov_results[i,] <- vr$vcov
    results$nll_optimum[i] <- -as.numeric(summary(lme4_model)$logLik)
    results$num_optim_iter[i] <- lme4_model@optinfo$feval
  })
}
log_real_world_result("lme4", results$time_estimation[i])

i <- i + 1
additional <- append_additional_results(results, cov_results, beta_results, i,
                                        Y, X, group_data, "gaussian",
                                        init_cov_pars, init_betas)
results <- additional$results
cov_results <- additional$cov_results
beta_results <- additional$beta_results
log_real_world_summary(results)

################################################################################
all_results <- cbind(results, cov_results, beta_results)
save_real_world_results(all_results, data_info, "./../results/chicago_building_permits.rds")
