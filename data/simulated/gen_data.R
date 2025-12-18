
make_data <- function(n, m1, sigma2=0.5^2, sigma2_1=0.5^2, sigma2_2=0.5^2,
                      likelihood, factor_m2=1, balanced = TRUE, size_neg_bin = 1,
                      has_F = FALSE, num_covariates=5, randef = NULL){
  if (!balanced) {
    if (size_neg_bin < 0 ) {
      group <- rpois(n = m1, lambda = n/m1 - 1) + 1
    } else {
      group <- rnbinom(n = m1, size = size_neg_bin, mu = n/m1 - 1) + 1
    }
    group <- rep(seq_along(group), times = group)
  } else {
    group <- rep(1,n)
    for(i in 1:m1) group[((i-1)*n/m1+1):(i*n/m1)] <- i
  }
  m2 <- as.integer(factor_m2 * m1)
  if (!balanced) {
    if (size_neg_bin < 0 ) {
      group2 <- rpois(n = m2, lambda = n/m2 - 1) + 1
    } else {
      group2 <- rnbinom(n = m2, size = size_neg_bin, mu = n/m2 - 1) + 1
    }
    group2 <- rep(seq_along(group2), times = group2)
  } else {
    group2 <- rep(1,n)
    for(i in 1:m2) group2[((i-1)*n/m2+1):(i*n/m2)] <- i
  }
  group2 <- group2[sample.int(n=length(group2), size=length(group2), replace=FALSE)]
  if (!balanced) {
    # make sure that both random effects have the same sample size by adding some samples
    if(length(group) < length(group2)) {
      group <- c(group, group[sample.int(n=length(group), size=length(group2) - length(group), replace=FALSE)])
    } else if(length(group) > length(group2)) {
      group2 <- c(group2, group2[sample.int(n=length(group2), size=length(group) - length(group2), replace=FALSE)])
    }
    n <- length(group)
  }
  b1 <- sqrt(sigma2_1) * rnorm(length(group))
  b2 <- sqrt(sigma2_2) * rnorm(length(unique(group2)))
  if (length(unique(group2)) != max(group2)) stop("not all levels samples -> gives index problems")
  eps <- b1[group] + b2[group2]
  group_data <- cbind(group,group2)
  #Simulate fixed effects
  if (has_F) {
    #The covariates X are sampled from a normal distribution with mean 0 and variance chosen such that 
    #the signal-to-noise ratio between the fixed and random effects (sigma2) is one, 
    #and the true regression coefficients are all 1 except for the intercept which is 0.
    beta <- c(0,rep(1,num_covariates))
    X <- cbind(rep(1,n),matrix(rnorm(num_covariates*n),nrow=n))
    f <- as.vector(X%*%beta)
    delta_var_f <- sqrt(sigma2/var(f))
    X[,-1] <- X[,-1] * delta_var_f
    f <- as.vector(X%*%beta)
    colnames(X) <- c("1",paste0("Cov_",1:num_covariates))
  }else{
    X <- matrix(rep(0,(num_covariates+1)*n),ncol=num_covariates+1)
    colnames(X) <- c("1",paste0("Cov_",1:num_covariates))
    f <- 0
  }
  #Simulate response variable
  if (likelihood == "bernoulli_probit") {
    probs <- pnorm(f+eps)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(f+eps)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(f+eps)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    mu <- exp(f+eps)
    y <- qgamma(runif(n), scale = mu, shape = 1)
  } else if (likelihood == "gaussian") {
    mu <- f + eps
    y <- sqrt(sigma2) * rnorm(n) + mu
  }
  list(y=y, X=X, group_data=group_data, eps=eps, f=f)
}
