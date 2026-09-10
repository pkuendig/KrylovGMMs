# Real-world applications

To replicate the real-world data applications from Section 6, proceed as follows:

1. Extract the ZIP files in `.\..\data\real_world`. 

2. Run the following scripts to perform estimation on the different real-world data sets:
    - `amazon_employee_access.R`
    - `cars.R`
    - `chicago_building_permits.R`
    - `instEval.R`
    - `KDDCup09_upselling.R`
    - `MovieLens.R`

3. Run `show_results.R` to generate the tables and plots from Section 6.

The application scripts also fit MixedModels.jl and R-INLA. The R-INLA fit uses
the empirical-Bayes integration strategy (`int.strategy = "eb"`) with flat
fixed-effect and log-precision priors. Since R-INLA integrates the fixed effects
out as part of the latent field, its empirical-Bayes optimum is the REML and not
the ML optimum: the reported variance parameters correspond to lme4 with
`REML = TRUE`, and the reported `nll_optimum` is a REML-type criterion that is
not on the same scale as the values of the other methods. The latent-field
marginals are not used and are switched off so that the measured runtime
reflects parameter estimation only. MixedModels.jl uses maximum likelihood
for Gaussian responses (`REML = false`) and a one-point Laplace approximation
for Bernoulli responses, matching the lme4 comparison.

To select the methods to run, set the comma-separated environment variable
`KRYLOVGMM_METHODS` before launching a script. For example, the following trial
omits lme4 and glmmTMB:

```text
KRYLOVGMM_METHODS=cholesky,krylov,mixedmodels,inla Rscript cars.R
```

The equivalent R-session configuration is

```r
options(krylovgmms.methods = c("cholesky", "krylov", "mixedmodels", "inla"))
source("cars.R")
```

Valid identifiers are `cholesky`, `krylov`, `lme4`, `glmmtmb`, `mixedmodels`,
and `inla`. MixedModels.jl requires Julia, the Julia package `MixedModels`, and
the R package `JuliaCall`. R-INLA requires the `INLA` R package.

MixedModels.jl performs a small untimed warm-up fit by default, excluding Julia
JIT compilation from `time_estimation`. To include compilation time, set
`options(krylovgmms.mixedmodels.warmup = FALSE)` before sourcing an application
script.
