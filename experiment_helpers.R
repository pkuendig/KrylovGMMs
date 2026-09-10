# Shared configuration and adapters for the experiment scripts.

set_experiment_working_directory <- function() {
  path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                   error = function(e) "")
  if (!nzchar(path)) {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) path <- sub("^--file=", "", file_arg[1])
  }
  if (nzchar(path)) setwd(dirname(normalizePath(path)))
  invisible(path)
}

experiment_method_ids <- c(
  "cholesky", "krylov", "lme4", "glmmtmb", "mixedmodels", "inla"
)

get_experiment_methods <- function() {
  methods <- getOption("krylovgmms.methods", NULL)
  if (is.null(methods)) {
    value <- Sys.getenv("KRYLOVGMM_METHODS", unset = "")
    methods <- if (nzchar(value)) strsplit(value, ",", fixed = TRUE)[[1]] else experiment_method_ids
  }
  methods <- tolower(trimws(methods))
  aliases <- c(
    gpboost_cholesky = "cholesky", gpboost_krylov = "krylov",
    mixedmodelsjl = "mixedmodels", `r-inla` = "inla"
  )
  methods <- ifelse(methods %in% names(aliases), aliases[methods], methods)
  unknown <- setdiff(methods, experiment_method_ids)
  if (length(unknown)) {
    stop("Unknown experiment method(s): ", paste(unknown, collapse = ", "),
         ". Valid values are: ", paste(experiment_method_ids, collapse = ", "))
  }
  unique(methods)
}

experiment_methods <- get_experiment_methods()
method_enabled <- function(method) method %in% experiment_methods

experiment_method_labels <- c(
  cholesky = "Cholesky (GPBoost)", krylov = "Krylov (GPBoost)",
  lme4 = "lme4", glmmtmb = "glmmTMB", mixedmodels = "MixedModels.jl",
  inla = "R-INLA (EB)"
)

experiment_log <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

log_time_vs_m_start <- function(script_name, likelihood, design) {
  labels <- unname(experiment_method_labels[experiment_methods])
  experiment_log("Starting ", script_name, ": ", likelihood,
                 " likelihood, ", design, " random-effects design; selected methods: ",
                 paste(labels, collapse = ", "))
}

log_time_vs_m_method <- function(m, method) {
  experiment_log("m=", format(m, big.mark = ",", scientific = FALSE),
                 ": running ", experiment_method_labels[[method]])
}

log_time_vs_m_iteration <- function(m, likelihood, design) {
  message("")
  experiment_log("Experiment: ", likelihood, " likelihood, ", design,
                 " random-effects design; m=",
                 format(m, big.mark = ",", scientific = FALSE))
}

log_time_vs_m_summary <- function(m, iteration_results) {
  label_to_id <- setNames(names(experiment_method_labels), experiment_method_labels)
  status <- vapply(seq_len(nrow(iteration_results)), function(j) {
    label <- iteration_results$method[j]
    method <- unname(label_to_id[label])
    elapsed <- iteration_results$time_estimation[j]
    state <- if (!method %in% experiment_methods) {
      "disabled"
    } else if (is.na(elapsed)) {
      "not run (size cutoff or failure)"
    } else if (!is.finite(elapsed)) {
      "failed or timed out"
    } else {
      paste0(format(round(elapsed, 3), nsmall = 3), " s")
    }
    paste0(label, " = ", state)
  }, character(1))
  experiment_log("m=", format(m, big.mark = ",", scientific = FALSE),
                 " complete: ", paste(status, collapse = "; "))
}

.current_real_world_dataset <- NA_character_

current_real_world_dataset <- function() {
  if (is.na(.current_real_world_dataset) || !nzchar(.current_real_world_dataset)) {
    "unknown dataset"
  } else {
    .current_real_world_dataset
  }
}

log_real_world_start <- function(dataset_name) {
  .current_real_world_dataset <<- dataset_name
  labels <- unname(experiment_method_labels[experiment_methods])
  experiment_log("Starting real-world experiment for dataset '", dataset_name, "'",
                 "; selected methods: ", paste(labels, collapse = ", "))
}

log_real_world_method <- function(method) {
  experiment_log("Dataset '", current_real_world_dataset(), "': running ",
                 experiment_method_labels[[method]])
}

log_real_world_result <- function(method, elapsed) {
  state <- if (!method %in% experiment_methods) {
    "disabled"
  } else if (is.na(elapsed)) {
    "failed or not run"
  } else if (!is.finite(elapsed)) {
    "failed or timed out"
  } else {
    paste0("completed in ", format(round(elapsed, 3), nsmall = 3), " s")
  }
  experiment_log("Dataset '", current_real_world_dataset(), "': ",
                 experiment_method_labels[[method]], " = ", state)
}

log_real_world_summary <- function(results) {
  label_to_id <- setNames(names(experiment_method_labels), experiment_method_labels)
  status <- vapply(seq_len(nrow(results)), function(j) {
    label <- results$method[j]
    method <- unname(label_to_id[label])
    elapsed <- results$time_estimation[j]
    state <- if (!method %in% experiment_methods) {
      "disabled"
    } else if (is.na(elapsed)) {
      "failed or not run"
    } else if (!is.finite(elapsed)) {
      "failed or timed out"
    } else {
      paste0(format(round(elapsed, 3), nsmall = 3), " s")
    }
    paste0(label, " = ", state)
  }, character(1))
  experiment_log("Dataset '", current_real_world_dataset(),
                 "' complete: ", paste(status, collapse = "; "))
}

# Methods that were not run in this session (i.e. that are not selected via
# `krylovgmms.methods` / KRYLOVGMM_METHODS) keep the values that are already
# stored in the results file, so that running a subset of the methods does not
# discard the results of the others. Rows are matched by method label.
merge_stored_real_world_results <- function(all_results, file) {
  if (!file.exists(file)) return(all_results)
  previous <- tryCatch(readRDS(file)$all_results, error = function(e) NULL)
  if (is.null(previous) || !"method" %in% names(previous)) return(all_results)

  label_to_id <- setNames(names(experiment_method_labels), experiment_method_labels)
  shared_cols <- intersect(names(all_results), names(previous))
  dropped_cols <- setdiff(names(previous), names(all_results))
  kept <- character(0)

  for (row in seq_len(nrow(all_results))) {
    label <- all_results$method[row]
    method <- unname(label_to_id[label])
    if (is.na(method) || method_enabled(method)) next
    old_row <- which(previous$method == label)
    if (length(old_row) != 1) next
    all_results[row, shared_cols] <- previous[old_row, shared_cols]
    all_results[row, setdiff(names(all_results), shared_cols)] <- NA
    kept <- c(kept, label)
  }

  if (length(kept)) {
    experiment_log("Dataset '", current_real_world_dataset(),
                   "': keeping stored results for ", paste(kept, collapse = ", "))
    if (length(dropped_cols)) {
      warning("Stored results contain columns that the current run does not: ",
              paste(dropped_cols, collapse = ", "),
              ". They are dropped for the methods that were not re-run.",
              call. = FALSE)
    }
  }
  all_results
}

# A missing or non-finite log-likelihood at the optimum means that the method did
# not produce a usable fit (it crashed, returned NaNs, or stopped before
# converging). Such a row is blanked out completely so that no runtime and no
# parameter estimates are reported for it.
invalidate_incomplete_results <- function(all_results) {
  if (!"nll_optimum" %in% names(all_results)) return(all_results)
  incomplete <- !is.finite(all_results$nll_optimum) &
    !apply(all_results[, setdiff(names(all_results), "method"), drop = FALSE], 1,
           function(row) all(is.na(row)))
  if (any(incomplete)) {
    experiment_log("Dataset '", current_real_world_dataset(),
                   "': no log-likelihood at the optimum for ",
                   paste(all_results$method[incomplete], collapse = ", "),
                   " -- runtime and estimates are set to NA")
    all_results[incomplete, setdiff(names(all_results), "method")] <- NA
  }
  all_results
}

save_real_world_results <- function(all_results, data_info, file) {
  all_results <- merge_stored_real_world_results(all_results, file)
  all_results <- invalidate_incomplete_results(all_results)
  saveRDS(list(all_results = all_results, data_info = data_info), file)
  invisible(all_results)
}

# Number of nonzeros in Z'Z for random-intercept incidence matrices, without
# constructing Z or requiring lme4. Diagonal blocks contribute one entry per
# level; each off-diagonal block contributes one per observed level pair.
count_group_ztz_nnz <- function(group_data) {
  groups <- as.data.frame(group_data)
  result <- sum(vapply(groups, function(x) length(unique(x)), integer(1)))
  if (ncol(groups) > 1) {
    pairs <- combn(seq_len(ncol(groups)), 2)
    result <- result + 2 * sum(apply(pairs, 2, function(j) {
      nrow(unique(groups[, j, drop = FALSE]))
    }))
  }
  result
}

# lme4 and glmmTMB do not parameterize the random-effect variances the way
# GPBoost does, so the initial values have to be transformed for all methods to
# start from the same point. `init_cov_pars` is (error variance, random-effect
# variances) for Gaussian likelihoods and (random-effect variances) otherwise.
#
# lme4's theta is the relative covariance factor, i.e. SD_j / sigma for lmer and
# SD_j for glmer (where sigma is fixed to one). lme4 also sorts the
# random-effect terms by decreasing number of levels via rev(order(nlev)) --
# which reverses ties relative to the formula order -- so theta has to be
# permuted into that order. Note that only glmer accepts a `fixef` starting
# value; lmer profiles the fixed effects out of its objective and rejects it.
lme4_start_theta <- function(init_cov_pars, group_data, likelihood) {
  variances <- if (identical(likelihood, "gaussian")) {
    init_cov_pars[-1] / init_cov_pars[1]
  } else {
    init_cov_pars
  }
  nlev <- vapply(as.data.frame(group_data), function(x) length(unique(x)), numeric(1))
  unname(sqrt(variances)[rev(order(nlev))])
}

# glmmTMB's theta is log(SD_j), in formula order, and the Gaussian residual SD is
# the separate parameter betadisp = log(sigma).
glmmtmb_start_theta <- function(init_cov_pars, likelihood) {
  variances <- if (identical(likelihood, "gaussian")) init_cov_pars[-1] else init_cov_pars
  unname(log(sqrt(variances)))
}

glmmtmb_start_betadisp <- function(init_cov_pars) unname(log(sqrt(init_cov_pars[1])))

require_experiment_package <- function(package, method = package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Method '", method, "' requires the R package '", package, "'.", call. = FALSE)
  }
}

.mixedmodels_initialized <- FALSE
.mixedmodels_warmed <- new.env(parent = emptyenv())

initialize_mixedmodels <- function() {
  if (.mixedmodels_initialized) return(invisible(NULL))
  require_experiment_package("JuliaCall", "mixedmodels")
  JuliaCall::julia_setup(installJulia = FALSE)
  JuliaCall::julia_library("MixedModels")
  JuliaCall::julia_library("LinearAlgebra")
  JuliaCall::julia_eval("
    begin
    function krylovgmms_fit(y, X, group_values, fixed_names, group_names, likelihood,
                            init_cov_pars, init_betas)
      # JuliaCall simplifies length-one R vectors to scalars. Normalize them
      # here so one-column designs and one grouping factor follow the same
      # code path as their multi-column counterparts.
      fixed_names = fixed_names isa AbstractString ? [fixed_names] : fixed_names
      group_names = group_names isa AbstractString ? [group_names] : group_names
      init_cov_pars = init_cov_pars isa Number ? [init_cov_pars] : init_cov_pars
      init_betas = init_betas isa Number ? [init_betas] : init_betas

      column_names = Symbol[:y]
      columns = Any[collect(y)]
      for (j, nm) in enumerate(fixed_names)
        push!(column_names, Symbol(nm))
        push!(columns, X[:, j])
      end
      for (j, nm) in enumerate(group_names)
        push!(column_names, Symbol(nm))
        push!(columns, string.(group_values[:, j]))
      end
      data = NamedTuple{Tuple(column_names)}(Tuple(columns))

      rhs = MixedModels.StatsModels.term(0)
      for nm in fixed_names
        rhs = rhs + MixedModels.StatsModels.term(Symbol(nm))
      end
      for nm in group_names
        rhs = rhs + (MixedModels.StatsModels.term(1) |
                     MixedModels.StatsModels.term(Symbol(nm)))
      end
      frm = MixedModels.StatsModels.term(:y) ~ rhs

      if likelihood == \"gaussian\"
        model = LinearMixedModel(frm, data)
        # MixedModels parameterizes random-effect SDs relative to the residual SD.
        theta0 = sqrt.(init_cov_pars[2:end] ./ init_cov_pars[1])
        model.optsum.initial .= krylovgmms_reorder(model, group_names, theta0)
        fit!(model; REML=false, progress=false)
      else
        model = GeneralizedLinearMixedModel(frm, data, Bernoulli())
        theta0 = krylovgmms_reorder(model, group_names, sqrt.(init_cov_pars))
        model.beta .= init_betas
        model.theta .= theta0
        # fit!(...; fast=false) builds its starting point as vcat(beta, optsum.final).
        model.optsum.final .= theta0
        fit!(model; nAGQ=1, progress=false)
      end
      return model
    end

    # MixedModels sorts the random-effect terms by number of levels (descending),
    # so a vector given in group_names order has to be permuted before it is
    # written to optsum.initial / model.theta. Valid because every term here is a
    # scalar random intercept, i.e. exactly one theta per grouping factor.
    function krylovgmms_reorder(model, group_names, values)
      names_in_model = String.(MixedModels.fname.(model.reterms))
      by_name = Dict(zip(group_names, values))
      return [by_name[name] for name in names_in_model]
    end

    function krylovgmms_re_variances(model, group_names, gaussian)
      group_names = group_names isa AbstractString ? [group_names] : group_names
      scale2 = gaussian ? MixedModels.varest(model) : 1.0
      names_in_model = String.(MixedModels.fname.(model.reterms))
      values = [only(LinearAlgebra.diag(Matrix(lambda)))^2 * scale2 for lambda in model.lambda]
      by_name = Dict(zip(names_in_model, values))
      return [by_name[name] for name in group_names]
    end
    end
  ")
  .mixedmodels_initialized <<- TRUE
  invisible(NULL)
}

prepare_mixedmodels_data <- function(y, X, group_data) {
  data <- data.frame(y = as.numeric(y), X, check.names = FALSE)
  groups <- as.data.frame(group_data, check.names = FALSE)
  groups[] <- lapply(groups, as.character)
  data.frame(data, groups, check.names = FALSE)
}

fit_mixedmodels <- function(y, X, group_data, likelihood, init_cov_pars, init_betas) {
  initialize_mixedmodels()
  JuliaCall::julia_assign("kg_fixed_names", colnames(X))
  JuliaCall::julia_assign("kg_group_names", colnames(group_data))
  JuliaCall::julia_assign("kg_likelihood", likelihood)
  JuliaCall::julia_assign("kg_init_cov_pars", as.numeric(init_cov_pars))
  JuliaCall::julia_assign("kg_init_betas", as.numeric(init_betas))

  warmup_signature <- paste(likelihood, paste(colnames(X), collapse = "|"),
                            paste(colnames(group_data), collapse = "|"), sep = "::")
  do_warmup <- isTRUE(getOption("krylovgmms.mixedmodels.warmup", TRUE)) &&
    !exists(warmup_signature, envir = .mixedmodels_warmed, inherits = FALSE)
  if (do_warmup) {
    n_warm <- 200L
    row_id <- seq_len(n_warm)
    warm_X <- vapply(seq_len(ncol(X)), function(j) {
      sin(row_id * (j + 0.37)) + cos(row_id / (j + 1.19))
    }, numeric(n_warm))
    colnames(warm_X) <- colnames(X)
    intercept <- tolower(colnames(X)) %in% c("intercept", "(intercept)")
    warm_X[, intercept] <- 1
    warm_groups <- vapply(seq_len(ncol(group_data)), function(j) {
      (row_id - 1L) %% (7L + j) + 1L
    }, integer(n_warm))
    if (ncol(group_data) == 1L) warm_groups <- matrix(warm_groups, ncol = 1L)
    colnames(warm_groups) <- colnames(group_data)
    warm_y <- if (identical(likelihood, "gaussian")) {
      as.numeric(warm_X %*% seq(0.05, 0.05 * ncol(X), by = 0.05) +
                   sin(row_id * 0.73))
    } else {
      as.integer(row_id %% 2L)
    }

    experiment_log("Warming up MixedModels.jl for ", likelihood,
                   " (Julia compilation is excluded from timings)")
    JuliaCall::julia_assign("kg_warm_y", warm_y)
    JuliaCall::julia_assign("kg_warm_X", unname(warm_X))
    JuliaCall::julia_assign("kg_warm_group_values", unname(warm_groups))
    invisible(JuliaCall::julia_eval(
      "krylovgmms_fit(kg_warm_y, kg_warm_X, kg_warm_group_values,
                       kg_fixed_names, kg_group_names, kg_likelihood,
                       kg_init_cov_pars, kg_init_betas)"
    ))
    assign(warmup_signature, TRUE, envir = .mixedmodels_warmed)
  }

  JuliaCall::julia_assign("kg_y", as.numeric(y))
  JuliaCall::julia_assign("kg_X", unname(as.matrix(X)))
  JuliaCall::julia_assign("kg_group_values", unname(as.matrix(group_data)))

  elapsed <- system.time(JuliaCall::julia_eval(
    "kg_model = krylovgmms_fit(kg_y, kg_X, kg_group_values,
                                kg_fixed_names, kg_group_names,
                                kg_likelihood, kg_init_cov_pars, kg_init_betas)"
  ))[["elapsed"]]
  gaussian <- identical(likelihood, "gaussian")
  JuliaCall::julia_assign("kg_gaussian", gaussian)
  list(
    method = "MixedModels.jl",
    time_estimation = elapsed,
    num_optim_iter = as.numeric(JuliaCall::julia_eval("kg_model.optsum.feval")),
    nll_optimum = as.numeric(JuliaCall::julia_eval("MixedModels.objective(kg_model) / 2")),
    beta = as.numeric(JuliaCall::julia_eval("collect(coef(kg_model))")),
    cov_pars = c(
      if (gaussian) as.numeric(JuliaCall::julia_eval("MixedModels.varest(kg_model)")),
      as.numeric(JuliaCall::julia_eval(
        "krylovgmms_re_variances(kg_model, kg_group_names, kg_gaussian)"
      ))
    )
  )
}

fit_inla_eb <- function(y, X, group_data, likelihood, init_cov_pars, init_betas) {
  require_experiment_package("INLA", "inla")
  fixed_names <- colnames(X)
  group_names <- colnames(group_data)
  data <- data.frame(y = as.numeric(y), X, check.names = FALSE)
  groups <- as.data.frame(group_data, check.names = FALSE)
  groups[] <- lapply(groups, function(x) as.integer(factor(x)))
  data <- data.frame(data, groups, check.names = FALSE)

  re_initial <- if (identical(likelihood, "gaussian")) init_cov_pars[-1] else init_cov_pars
  flat_prior <- "expression: return(0);"
  re_terms <- vapply(seq_along(group_names), function(j) {
    sprintf(
      "f(%s, model='iid', hyper=list(prec=list(prior='%s', initial=%.17g)))",
      group_names[j], flat_prior, -log(re_initial[j])
    )
  }, character(1))
  formula <- as.formula(paste(
    "y ~ -1 +", paste(c(fixed_names, re_terms), collapse = " + ")
  ))
  control_family <- if (identical(likelihood, "gaussian")) {
    list(hyper = list(prec = list(prior = flat_prior,
                                  initial = -log(init_cov_pars[1]))))
  } else {
    list()
  }

  elapsed <- system.time(model <- INLA::inla(
    formula,
    family = if (identical(likelihood, "gaussian")) "gaussian" else "binomial",
    data = data,
    # Zero precision gives the fixed effects a flat prior. Together with flat
    # priors on log-precisions, the EB mode maximizes the (INLA-approximated)
    # marginal likelihood. Note that INLA integrates the fixed effects out as
    # part of the latent field, so this targets the REML rather than the ML
    # optimum; the difference is O(p/n) and the estimates below match lme4 with
    # REML = TRUE, not the ML fits of the other methods.
    control.fixed = list(mean = 0, prec = 0),
    control.family = control_family,
    # The variance parameters are estimated from the Laplace approximation of
    # pi(theta|y) irrespective of `strategy`, which only governs the marginals of
    # the latent field. Those marginals are not used here, so the cheapest
    # strategy is used and their computation is switched off; this was verified
    # to leave the hyperparameter mode and mlik unchanged.
    control.inla = list(strategy = "gaussian", int.strategy = "eb"),
    control.compute = list(mlik = TRUE, return.marginals = FALSE),
    control.predictor = list(compute = FALSE),
    verbose = FALSE
  ))[["elapsed"]]

  # The EB point estimate is the mode of the internal (log-precision)
  # parameterization found by INLA's optimizer. The "mode" column of
  # summary.hyperpar is the mode of the marginal density after the nonlinear
  # transformation to the precision scale, which differs substantially.
  theta_mode <- model$mode$theta
  re_names <- paste("Log precision for", group_names)
  cov_pars <- exp(-as.numeric(theta_mode[re_names]))
  if (identical(likelihood, "gaussian")) {
    cov_pars <- c(exp(-as.numeric(
      theta_mode["Log precision for the Gaussian observations"])), cov_pars)
  }
  if (anyNA(cov_pars)) {
    stop("Could not match INLA hyperparameters: ",
         paste(names(theta_mode), collapse = ", "), call. = FALSE)
  }
  list(
    method = "R-INLA (EB)",
    time_estimation = elapsed,
    num_optim_iter = NA_real_,
    # Negative log marginal likelihood with the fixed effects integrated out, so
    # this is a REML-type criterion and is not on the same scale as the ML
    # values reported by the other methods.
    nll_optimum = -as.numeric(model$mlik[1, 1]),
    beta = as.numeric(model$summary.fixed[fixed_names, "mode"]),
    cov_pars = cov_pars
  )
}

fit_additional_method <- function(method, y, X, group_data, likelihood,
                                  init_cov_pars, init_betas) {
  tryCatch(
    switch(
      method,
      mixedmodels = fit_mixedmodels(y, X, group_data, likelihood, init_cov_pars, init_betas),
      inla = fit_inla_eb(y, X, group_data, likelihood, init_cov_pars, init_betas),
      stop("Unsupported additional method: ", method)
    ),
    error = function(e) {
      experiment_log("Dataset '", current_real_world_dataset(), "': method ",
                     experiment_method_labels[[method]], " FAILED: ",
                     conditionMessage(e))
      warning("Method '", method, "' failed: ", conditionMessage(e), call. = FALSE)
      list(
        method = c(mixedmodels = "MixedModels.jl", inla = "R-INLA (EB)")[[method]],
        time_estimation = NA_real_, num_optim_iter = NA_real_, nll_optimum = NA_real_,
        beta = rep(NA_real_, ncol(X)),
        cov_pars = rep(NA_real_, ncol(group_data) + as.integer(likelihood == "gaussian"))
      )
    }
  )
}

store_additional_results <- function(results, cov_results, beta_results, i, fit) {
  results$method[i] <- fit$method
  results$time_estimation[i] <- fit$time_estimation
  results$num_optim_iter[i] <- fit$num_optim_iter
  results$nll_optimum[i] <- fit$nll_optimum
  cov_results[i, ] <- fit$cov_pars
  beta_results[i, ] <- fit$beta
  list(results = results, cov_results = cov_results, beta_results = beta_results)
}

append_additional_results <- function(results, cov_results, beta_results, i,
                                      y, X, group_data, likelihood,
                                      init_cov_pars, init_betas) {
  labels <- c(mixedmodels = "MixedModels.jl", inla = "R-INLA (EB)")
  for (method in names(labels)) {
    if (method_enabled(method)) {
      log_real_world_method(method)
      fit <- fit_additional_method(method, y, X, group_data, likelihood,
                                   init_cov_pars, init_betas)
    } else {
      fit <- list(
        method = labels[[method]], time_estimation = NA_real_,
        num_optim_iter = NA_real_, nll_optimum = NA_real_,
        beta = rep(NA_real_, ncol(X)),
        cov_pars = rep(NA_real_, ncol(group_data) + as.integer(likelihood == "gaussian"))
      )
    }
    stored <- store_additional_results(results, cov_results, beta_results, i, fit)
    results <- stored$results
    cov_results <- stored$cov_results
    beta_results <- stored$beta_results
    log_real_world_result(method, fit$time_estimation)
    i <- i + 1
  }
  list(results = results, cov_results = cov_results, beta_results = beta_results, i = i)
}
