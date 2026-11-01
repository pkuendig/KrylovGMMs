# Simulation studies

Simulated data is generated using the function `make_data()` in the file `root\data\simulated\gen_data.R`.

* `bias_analysis.R`
  
   Variance parameter estimation with a Laplace approximation and Krylov subspace methods for varying numbers of repeated observations per random effect realization. This additional analysis is not reported in the article.

* `compare_preconditioner_nll.R`

  Preconditioner comparison with regard to runtime and accuracy of log-marginal likelihoods (Figure 1 in Section 5.2 and Figures A1-A5 in Appendix H).

* `estimation_prediction_Gaussian.R`, `estimation_prediction_Bernoulli.R`, and `estimation_prediction_show_results.R`
  
  Parameter estimation and prediction on 100 simulated data sets for Gaussian and Bernoulli likelihoods using Krylov (GPBoost), Cholesky (GPBoost), lme4, and glmmTMB. First run `estimation_prediction_Gaussian.R` and `estimation_prediction_Bernoulli.R` and then run `estimation_prediction_show_results.R` to generate the results of Section 5.5 (Figures A7 and A8 and Tables A1-A4 in Appendix J).

* `eval_predictive_variances.R`
  
  Compare simulation-based predictive variances from Algorithm 1 with Cholesky-based computations in terms of accuracy, measured by RMSE, and runtime (Figure 2 in Section 5.3). Results for the other approaches cannot be replicated as these are not implemented in the publicly available version of GPBoost.

* `time_vs_m_balanced_Gaussian.R`, `time_vs_m_balanced_Bernoulli.R`, and `time_vs_m_balanced_plotting.R`

  Runtime comparison for parameter estimation for different numbers of random effects m, when the random effects design is balanced (m_1=m_2). First run `time_vs_m_balanced_Gaussian.R` and `time_vs_m_balanced_Bernoulli.R` and then run `time_vs_m_balanced_plotting.R` to generate Figure 3 in Section 5.4.

* `time_vs_m_unbalanced_Gaussian.R`, `time_vs_m_unbalanced_Bernoulli.R`, and `time_vs_m_unbalanced_plotting.R`

  Runtime comparison for parameter estimation for different numbers of random effects m, when the random effects design is unbalanced (m_2=m_1/2). First run `time_vs_m_unbalanced_Gaussian.R` and `time_vs_m_unbalanced_Bernoulli.R` and then run `time_vs_m_unbalanced_plotting.R` to generate Figure A6 in Appendix J.

The `time_vs_m` experiments include MixedModels.jl and R-INLA in addition to
GPBoost, lme4, and glmmTMB. R-INLA uses its empirical-Bayes hyperparameter
integration strategy (`int.strategy = "eb"`) with flat fixed-effect and
log-precision priors. Because R-INLA integrates the fixed effects out as part of
the latent field, its empirical-Bayes optimum is the REML rather than the ML
optimum; the latent-field marginals are not used and are switched off so that
the measured runtime reflects parameter estimation only.

Methods can be selected without editing the scripts. For example, a trial run
that excludes lme4 and glmmTMB can be launched from a shell with

```text
KRYLOVGMM_METHODS=cholesky,krylov,mixedmodels,inla Rscript time_vs_m_balanced_Gaussian.R
```

or from an R session with

```r
options(krylovgmms.methods = c("cholesky", "krylov", "mixedmodels", "inla"))
source("time_vs_m_balanced_Gaussian.R")
```

Valid method identifiers are `cholesky`, `krylov`, `lme4`, `glmmtmb`,
`mixedmodels`, and `inla`.

By default, the MixedModels.jl adapter performs one small untimed warm-up fit
per likelihood/formula shape so that Julia JIT compilation is not included in
`time_estimation`. Set `options(krylovgmms.mixedmodels.warmup = FALSE)` before
running a script if startup and compilation time should be included instead.
