################################################################################
# Plot runtime results for parameter estimation for a balanced random effects design.
################################################################################

library(ggplot2)
library(ggh4x)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

Gaussian_results <- readRDS("./../results/time_vs_m_balanced_Gaussian.rds")
Bernoulli_results <- readRDS("./../results/time_vs_m_balanced_Bernoulli.rds")

Gaussian_results$likelihood <- "Gaussian"
Bernoulli_results$likelihood <- "Bernoulli"
results <- rbind(Gaussian_results, Bernoulli_results)

#Remove not evaluated combinations
results <- results[!is.na(results$time_estimation),]

method_levels <- c("Krylov (GPBoost)", "Cholesky (GPBoost)", "lme4", "glmmTMB",
                   "MixedModels.jl", "R-INLA (EB)")
results$method <- factor(results$method,
                         levels = intersect(method_levels, unique(results$method)))

# Fixed colors and shapes per method so that all figures below are comparable
# (first four colors / shapes are the ggplot defaults of ColorBrewer "Set1").
method_colors <- c("Krylov (GPBoost)"   = "#E41A1C",
                   "Cholesky (GPBoost)" = "#377EB8",
                   "lme4"               = "#4DAF4A",
                   "glmmTMB"            = "#984EA3",
                   "MixedModels.jl"     = "#FF7F00",
                   "R-INLA (EB)"        = "#A65628")
method_shapes <- c("Krylov (GPBoost)"   = 16,
                   "Cholesky (GPBoost)" = 17,
                   "lme4"               = 15,
                   "glmmTMB"            = 3,
                   "MixedModels.jl"     = 7,
                   "R-INLA (EB)"        = 8)
results$likelihood <- factor(results$likelihood, levels = c("Gaussian", "Bernoulli"))

# Axis ticks used in the runtime figure.
x_breaks <- c(1000, 2000, 3000, 5000, 10000, 20000,
              50000, 100000, 200000, 500000, 1000000)
bernoulli_y_breaks <- gaussian_y_breaks <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 4000, 8000)

x_labels <- function(x) format(x, big.mark = "'", scientific = FALSE, trim = TRUE)
gaussian_y_labels <- function(x) formatC(x, format = "f", digits = 1, big.mark = "'")
bernoulli_y_labels <- function(x) formatC(x, format = "f", digits = 0, big.mark = "'")

make_plot <- function(methods) {
  dat <- droplevels(results[results$method %in% methods,])
  ggplot(data=dat, aes(x=m, y=time_estimation, color=method, shape=method)) +
    geom_line(linewidth=1) + geom_point(size=2) +
    facet_nested_wrap(~likelihood, nrow=1, drop=T, scales = "free") +
    xlab("m") + ylab("Time (s)") + theme_bw() +
    theme(legend.position = "top", legend.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) +
    facetted_pos_scales(
      y = list(
        likelihood == "Gaussian" ~ scale_y_continuous(
          trans = "log1p", breaks = gaussian_y_breaks, labels = gaussian_y_labels),
        likelihood == "Bernoulli" ~ scale_y_continuous(
          trans = "log1p", breaks = bernoulli_y_breaks, labels = bernoulli_y_labels)
      ),
      x = list(
        likelihood == "Gaussian" ~ scale_x_log10(breaks = x_breaks, labels = x_labels),
        likelihood == "Bernoulli" ~ scale_x_log10(breaks = x_breaks, labels = x_labels)
      )
    )
}

save_plot <- function(p, file_name) {
  ggsave(p, file = paste0("./../plots/", file_name, ".png"),
         width = 8.5, height = 5, dpi = 300)
}

# All methods
p <- make_plot(levels(results$method))
p
save_plot(p, "time_vs_m_balanced")

# lme4 and glmmTMB only
p_lme4_glmmTMB <- make_plot(c("lme4", "glmmTMB"))
p_lme4_glmmTMB
save_plot(p_lme4_glmmTMB, "time_vs_m_balanced_lme4_glmmTMB")

p_gpboost <- make_plot(c("lme4", "glmmTMB", "Krylov (GPBoost)", "Cholesky (GPBoost)"))
p_gpboost
save_plot(p_gpboost, "time_vs_m_balanced_lme4_glmmTMB_gpboost")

p_mixedmodels <- make_plot(c("lme4", "glmmTMB", "Krylov (GPBoost)", "Cholesky (GPBoost)", "MixedModels.jl"))
p_mixedmodels
save_plot(p_mixedmodels, "time_vs_m_balanced_lme4_glmmTMB_gpboost_mixedmodels")
