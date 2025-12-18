library(gpboost)
library(lme4)
library(glmmTMB)

run_with_timeout <- function(expr, seconds = 60*60*24*5) {
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  setTimeLimit(elapsed = seconds, transient = TRUE)
  force(expr)
}

data_info <- NULL

set.seed(1)

###Import data##################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
my_data <- read.csv("./../data/real_world/ml-32m/ratings.csv")

###Prepare data#################################################################

#Categorical variables

cat_vars <- c("userId", "movieId")
group_data <- my_data[,cat_vars]

#Predictor variables
pred_vars <- c("time", "time_sq")
time <- my_data[,"timestamp"]
time <- (time - min(time)) / (max(time) - min(time))
X <- cbind(rep(1, dim(my_data)[1]),time, (time-mean(time))^2)
colnames(X) <- c("Intercept", "time", "time_sq")

#Response
Y <- my_data$rating
# summary(Y)
rm(my_data)
gc()

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

# Load previous results
# res_loaded <- readRDS("./../results/MovieLens_32m.rds")
# results <- res_loaded$all_results[,1:4]
# cov_results <- res_loaded$all_results[,5:7]
# beta_results <- res_loaded$all_results[,8:10]

# We use the same initial values for all libraries.
# The following are internal default initial values used by GPBoost
init_cov_pars <- c(0.560726, rep(0.280363, ncol(group_data)))
init_betas <- c(3.5404, rep(0, ncol(X)-1))

###Krylov#######################################################################
i <- 1
it_model <- GPModel(group_data = group_data,
                    likelihood="gaussian",
                    matrix_inversion_method = "iterative")
it_model$set_optim_params(params = list(maxit=1000,
                                        trace=TRUE,
                                        init_cov_pars=init_cov_pars,
                                        init_coef=init_betas,
                                        cg_preconditioner_type="ssor"))
results$method[i] <- "Krylov (GPBoost)"
time_result <- tryCatch(
  {
    system.time({
      run_with_timeout({
        inner_fit <- tryCatch(
          it_model$fit(y = Y, X = X),
          error = function(e) {
            message("❌ Model fitting crashed: ", conditionMessage(e))
            return(NA)
          }
        )
        inner_fit
      })
    })[["elapsed"]]
  },
  error = function(e) {
    message("⏰ Timeout or outer error: ", conditionMessage(e))
    return(Inf)
  }
)
results$time_estimation[i] <- time_result
if (results$time_estimation[i] != Inf & !is.na(results$time_estimation[i])) {
  beta_results[i,] <- it_model$get_coef()
  cov_results[i,] <- it_model$get_cov_pars()
  results$nll_optimum[i] <- it_model$get_current_neg_log_likelihood()
  results$num_optim_iter[i] <- it_model$get_num_optim_iter()
}

i <- i + 1

all_results <- cbind(results, cov_results, beta_results)
saveRDS(list(all_results=all_results, data_info=data_info), "./../results/MovieLens_32m.rds")
rm(it_model)
gc()

###Cholesky#####################################################################
## Crashes
chol_model <- GPModel(group_data = group_data,
                      likelihood="gaussian",
                      matrix_inversion_method = "cholesky")
chol_model$set_optim_params(params = list(maxit=1000,
                                          trace=TRUE,
                                          init_cov_pars=init_cov_pars,
                                          init_coef=init_betas))
results$method[i] <- "Cholesky (GPBoost)"
time_result <- tryCatch(
  {
    system.time({
      run_with_timeout({
        inner_fit <- tryCatch(
          chol_model$fit(y = Y, X = X),
          error = function(e) {
            message("❌ Model fitting crashed: ", conditionMessage(e))
            return(NA)
          }
        )
        inner_fit
      })
    })[["elapsed"]]
  },
  error = function(e) {
    message("⏰ Timeout or outer error: ", conditionMessage(e))
    return(Inf)
  }
)
results$time_estimation[i] <- time_result
if (results$time_estimation[i] != Inf & !is.na(results$time_estimation[i])) {
  beta_results[i,] <- chol_model$get_coef()
  cov_results[i,] <- chol_model$get_cov_pars()
  results$nll_optimum[i] <- chol_model$get_current_neg_log_likelihood()
  results$num_optim_iter[i] <- chol_model$get_num_optim_iter()
}

i <- i + 1

all_results <- cbind(results, cov_results, beta_results)
saveRDS(list(all_results=all_results, data_info=data_info), "./../results/MovieLens_32m.rds")
rm(chol_model)
gc()

###glmmTMB######################################################################
## Crashes (memory error)
formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
results$method[i] <- "glmmTMB"
glmmTMB_model <- NA
time_result <- tryCatch(
  {
    system.time({
      run_with_timeout({
        inner_fit <- tryCatch(
          glmmTMB_model <<- glmmTMB(formula, family=gaussian(), data = data.frame(y = Y, cbind(X, group_data)), 
                                   start=list(theta = init_cov_pars[-1], beta = init_betas)),
          error = function(e) {
            message("❌ Model fitting crashed: ", conditionMessage(e))
            return(NA)
          }
        )
        inner_fit
      })
    })[["elapsed"]]
  },
  error = function(e) {
    message("⏰ Timeout or outer error: ", conditionMessage(e))
    return(Inf)
  }
)
results$time_estimation[i] <- time_result
if (results$time_estimation[i] != Inf & !is.na(results$time_estimation[i])) {
  beta_results[i,] <- fixef(glmmTMB_model)$cond
  cov_results[i,] <- c((summary(glmmTMB_model)$sigma)^2, as.numeric(VarCorr(glmmTMB_model)$cond))
  results$nll_optimum[i] <- -as.numeric(summary(glmmTMB_model)$logLik)
  results$num_optim_iter[i] <- glmmTMB_model$fit$iterations
}

i <- i + 1
all_results <- cbind(results, cov_results, beta_results)
saveRDS(list(all_results=all_results, data_info=data_info), "./../results/MovieLens_32m.rds")
rm(glmmTMB_model)
gc()

###lme4#########################################################################
formula <- as.formula(paste0("y ~ -1 + ",paste0(colnames(X), collapse = ' + ')," + ",paste0("(1|",cat_vars,")", collapse = ' + '))) 
results$method[i] <- "lme4"
lme4_model <- NA
time_result <- tryCatch(
  {
    system.time({
      run_with_timeout({
        inner_fit <- tryCatch(
          lme4_model <<- lmer(formula, REML = FALSE, data = data.frame(y = Y, cbind(X, group_data)), 
                             start=list(theta = init_cov_pars[-1], fixef = init_betas)),
          error = function(e) {
            message("❌ Model fitting crashed: ", conditionMessage(e))
            return(NA)
          }
        )
        inner_fit
      })
    })[["elapsed"]]
  },
  error = function(e) {
    message("⏰ Timeout or outer error: ", conditionMessage(e))
    return(Inf)
  }
)
results$time_estimation[i] <- time_result
if (results$time_estimation[i] != Inf & !is.na(results$time_estimation[i])) {
  beta_results[i,] <- summary(lme4_model)$coefficients[,1]
  vr <- as.data.frame(VarCorr(lme4_model))
  vr <- vr[match(c("Residual", cat_vars), vr$grp),]
  cov_results[i,] <- vr$vcov
  results$nll_optimum[i] <- -as.numeric(summary(lme4_model)$logLik)
  results$num_optim_iter[i] <- lme4_model@optinfo$feval
}

################################################################################
all_results <- cbind(results, cov_results, beta_results)
saveRDS(list(all_results=all_results, data_info=data_info), "./../results/MovieLens_32m.rds")


###Store data info##############################################################
data_cols <- c("ds_name", "n", "p", "K", "response_var", "class_imbalance", "cat_var", "nr_levels", "nnz_ZtZ")
data_info <- data.frame(matrix(nrow=length(cat_vars), ncol = length(data_cols)))
colnames(data_info) <- data_cols

data_info$ds_name[1] <- "MovieLens_32m"
data_info$n[1] <- length(Y)
data_info$p[1] <- ncol(X) - 1
data_info$K[1] <- ncol(group_data)
data_info$response_var[1] <- "rating"

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
#image(ZtZ, xlab="", ylab="", sub ="", main="MovieLens_32m")
rm(ZtZ)
rm(parsedFormula)
gc()
