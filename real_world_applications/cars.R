library(gpboost)
library(lme4)
library(glmmTMB)

set.seed(1)

###Import data##################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
my_data <- read.csv("./../data/real_world/cars_df5.csv")

###Prepare data#################################################################

#Categorical variables
cat_vars <- c("model_id", "location_id")
group_data_df <- my_data[,cat_vars]
group_data <- matrix(nrow = nrow(group_data_df), ncol=ncol(group_data_df))
colnames(group_data) <- colnames(group_data_df)
for(j in 1:dim(group_data)[2]) {
  group_data[,j] <- as.numeric(as.factor(group_data_df[,j]))
}

#Predictor variables
feat_cols <- colnames(my_data)
feat_cols <- feat_cols[!(feat_cols %in% c("price", cat_vars, 
                                          "manufacturerbmw", "conditionexcellent", "fueldiesel", "title_statusclean",
                                          "transmissionautomatic", "drive4wd", "sizecompact", "typebus", "paint_colorblack"))]
X <- as.matrix(my_data[,feat_cols])
X <- cbind(intercept = rep(1,nrow(X)), X)

#Response
Y <- log(my_data$price)

###Store data info#############################################################
data_cols <- c("ds_name", "n", "p", "K", "response_var", "class_imbalance", "cat_var", "nr_levels", "nnz_ZtZ")
data_info <- data.frame(matrix(nrow=length(cat_vars), ncol = length(data_cols)))
colnames(data_info) <- data_cols

data_info$ds_name[1] <- "cars"
data_info$n[1] <- length(Y)
data_info$p[1] <- ncol(X) - 1
data_info$K[1] <- ncol(group_data)
data_info$response_var[1] <- "log(price)"

for(i in 1:ncol(group_data)){
  data_info$cat_var[i] <- cat_vars[i]
  data_info$nr_levels[i] <- length(unique(group_data[,cat_vars[i]]))
}

#Infos on Z^TZ
parsedFormula <- lFormula(form = paste0("y ~ 1 + ", paste0("(1|",cat_vars,")", collapse = ' + ')), 
                          data = data.frame(cbind(y=rep(1,dim(group_data)[1]),group_data))) 
Ztlist <- parsedFormula$reTrms$Ztlist
for(i in 1:length(Ztlist)){
  if(i==1){
    Z <- t(Ztlist[[i]])
  }else{
    Z <- cbind(Z, t(Ztlist[[i]]))
  }
}
ZtZ <- t(Z)%*%Z
data_info$nnz_ZtZ[1] <- nnzero(ZtZ)
#image(ZtZ, xlab="", ylab="", sub="", main="cars")

###Data frames##################################################################
res_cols <- c("method", "time_estimation", "num_optim_iter", "nll_optimum")
results <- data.frame(matrix(nrow=4, ncol = length(res_cols)))
colnames(results) <- res_cols

cov_cols <- c("sigma", paste0("sigma2_", 1:length(cat_vars)))
cov_results <- data.frame(matrix(nrow=4, ncol = length(cov_cols)))
colnames(cov_results) <- cov_cols

beta_cols <- paste0("beta_", 0:(ncol(X)-1))
beta_results <- data.frame(matrix(nrow=4, ncol = length(beta_cols)))
colnames(beta_results) <- beta_cols

###Estimation###################################################################
i <- 1

# We use the same initial values for all libraries.
# The following are internal default initial values used by GPBoost
init_cov_pars <- c(0.282, 0.141, 0.141)
init_betas <- c(9.597, rep(0, ncol(X)-1))

###Cholesky#####################################################################
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

i <- i + 1

###Krylov#####################################################################
it_model <- GPModel(group_data = group_data,
                    likelihood="gaussian",
                    matrix_inversion_method = "iterative")

it_model$set_optim_params(params = list(maxit=1000,
                                        trace=TRUE,
                                        init_cov_pars=init_cov_pars,
                                        init_coef=init_betas,
                                        seed_rand_vec_trace=1,
                                        cg_preconditioner_type="zero_infill_incomplete_cholesky"))


results$method[i] <- "Krylov (GPBoost)"
results$time_estimation[i] <- system.time(it_model$fit(y=Y, X=X))[3]
beta_results[i,] <- it_model$get_coef()
cov_results[i,] <- it_model$get_cov_pars()
results$nll_optimum[i] <- it_model$get_current_neg_log_likelihood()
results$num_optim_iter[i] <- it_model$get_num_optim_iter()

i <- i + 1

###glmmTMB####################################################################
formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
results$time_estimation[i] <- system.time(glmmTMB_model <- glmmTMB(formula, family=gaussian(), data = data.frame(y = Y, cbind(X, group_data)), 
                                                                   start=list(theta = init_cov_pars[-1], beta = init_betas)))[3]

results$method[i] <- "glmmTMB"
beta_results[i,] <- fixef(glmmTMB_model)$cond
cov_results[i,] <- c((summary(glmmTMB_model)$sigma)^2, as.numeric(VarCorr(glmmTMB_model)$cond))
results$nll_optimum[i] <- -as.numeric(summary(glmmTMB_model)$logLik)
results$num_optim_iter[i] <- glmmTMB_model$fit$iterations

i <- i + 1

###lme4#########################################################################
try({
  formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
  results$time_estimation[i] <- system.time(lme4_model <- lmer(formula, REML = FALSE, data = data.frame(y = Y, cbind(X, group_data)), 
                                                               start=list(theta = init_cov_pars[-1], fixef = init_betas)))[3]

  results$method[i] <- "lme4"
  beta_results[i,] <- summary(lme4_model)$coefficients[,1]
  vr <- as.data.frame(VarCorr(lme4_model))
  vr <- vr[match(c("Residual", cat_vars), vr$grp),]
  cov_results[i,] <- vr$vcov
  results$nll_optimum[i] <- -as.numeric(summary(lme4_model)$logLik)
  results$num_optim_iter[i] <- lme4_model@optinfo$feval
})

################################################################################
all_results <- cbind(results, cov_results, beta_results)
saveRDS(list(all_results=all_results, data_info=data_info), "./cars.rds")
