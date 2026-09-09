# Scalable Krylov Subspace Methods for Generalized Mixed Effects Models with Crossed Random Effects

This repository contains the R code for the simulation studies and real-world applications in the paper "Scalable Krylov Subspace Methods for Generalized Mixed Effects Models with Crossed Random Effects". Krylov subspace methods for Generalized Mixed Effects Models (GMMs) are implemented in the package GPBoost: <https://github.com/fabsig/GPBoost>.

* The `simulation_studies` folder contains the scripts to run the simulations. See `simulation_studies\README.md` for a detailed description of the different simulations.
* The `real_world_applications` folder contains the scripts to run the real-world data applications. See `real_world_applications\README.md` for more information about the real-world applications.
* The `data` folder contains a file for generating simulated data and the data sets for the real-world applications. See `data\real_world\README.md` for more information on the sources of the real-world datasets.

### Structure of the repository
```
root
│   README.md
│
└───data
│   │
│   └───real_world
│   │   │   README.md
│   │   │   Amazon_employee_access.csv
│   │   │   cars_df5.csv
│   │   │   chicago_building_permits.zip
│   │   │   KDDCup09_upselling.zip
│   │
│   └───simulated
│       │   gen_data.R
│
└───real_world_applications
│   │   README.md
│   │   amazon_employee_access.R
│   │   cars.R
│   │   chicago_building_permits.R
│   │   instEval.R
│   │   KDDCup09_upselling.R
│   │   MovieLens.R
│   │   show_results.R
│
└───simulation_studies
    │   README.md
    │   bias_analysis.R
    │   compare_preconditioner_nll.R
    │   estimation_prediction_Bernoulli.R
    │   estimation_prediction_Gaussian.R
    │   estimation_prediction_show_results.R 
    │   eval_predictive_variances.R 
    │   time_vs_m_balanced_Bernoulli.R 
    │   time_vs_m_balanced_Gaussian.R
    │   time_vs_m_balanced_plotting.R
    │   time_vs_m_unbalanced_Bernoulli.R
    │   time_vs_m_unbalanced_Gaussian.R
    │   time_vs_m_unbalanced_plotting.R           
```
