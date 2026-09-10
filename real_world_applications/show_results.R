
library(xtable)
library(ggplot2)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_sets <- c("cars", "chicago_building_permits", "instEval", "MovieLens_32m",
               "amazon_employee_access", "KDDCup09_upselling")

# A missing or non-finite log-likelihood at the optimum means that the method did
# not produce a usable fit (it crashed, returned NaNs, or stopped before
# converging). Neither its runtime nor its estimates are reported. The
# experiment scripts already do this when saving (see experiment_helpers.R); it
# is repeated here so that older results files are treated the same way.
read_results <- function(d) {
  my_data <- readRDS(paste0("./../results/", d, ".rds"))
  incomplete <- !is.finite(my_data$all_results$nll_optimum)
  my_data$all_results[incomplete, setdiff(names(my_data$all_results), "method")] <- NA
  my_data
}

# A plain comma between digits is treated as mathematical punctuation in LaTeX
# and is typeset with a visible space (e.g. '100, 000'). '{,}' is therefore used
# as the thousands separator in all LaTeX output.
latex_format_args <- list(big.mark = "{,}")

###Table: data info#############################################################

#import results
info_data <- data.frame()
for(d in 1:length(data_sets)){
  my_data <- read_results(data_sets[d])
  if(d==1){
    info_data <- my_data$data_info
  }
  else{
    info_data <- rbind(info_data, my_data$data_info)
  }
}

print(xtable(info_data), format.args = latex_format_args)

###Table: runtimes##############################################################

#import results
results <- data.frame()
c <- c(1,2,4)
for(d in 1:length(data_sets)){
  my_data <- read_results(data_sets[d])
  new_results <-  my_data$all_results[,c]
  new_results <- cbind(ds_name = c(data_sets[d], rep(NA, nrow(new_results)-1)), new_results)
  if(d==1){
    results <- new_results
  }
  else{
    results <- rbind(results, new_results)
  }
}

print(xtable(results), format.args = latex_format_args)

###Plots########################################################################

#import results
results <- data.frame()
c <- c(1,2,4)
for(d in 1:length(data_sets)){
  my_data <- read_results(data_sets[d])
  new_results <-  my_data$all_results[,c]
  new_results$ds_name <- data_sets[d]
  if(d==1){
    results <- new_results
  }
  else{
    results <- rbind(results, new_results)
  }
}

# Methods with an NA runtime crashed (out of memory / no convergence). They are
# shown separately at the top of the runtime plot
crashed <- results[is.na(results$time_estimation),]
results <- results[!is.na(results$time_estimation),]

method_levels <- c("Krylov (GPBoost)", "Cholesky (GPBoost)", "lme4", "glmmTMB",
                   "MixedModels.jl", "R-INLA (EB)")
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

results$method <- factor(results$method,
                         levels = intersect(method_levels, unique(results$method)))
results$ds_name <- factor(results$ds_name, levels = data_sets)
crashed$method <- factor(crashed$method, levels = levels(results$method))
crashed$ds_name <- factor(crashed$ds_name, levels = data_sets)

# Shorter data set names for the x-axis
ds_labels <- c("cars"                     = "cars",
               "chicago_building_permits" = "building_permits",
               "instEval"                 = "instEval",
               "amazon_employee_access"   = "employee_access",
               "KDDCup09_upselling"       = "upselling",
               "MovieLens_32m"            = "MovieLens_32m")

y_labels <- function(x) format(x, big.mark = "'", scientific = FALSE, trim = TRUE)

save_plot <- function(p, file_name, width = 7, height = 7) {
  ggsave(p, file = paste0("./../plots/", file_name, ".png"),
         width = width, height = height, dpi = 300, device = ragg::agg_png)
}

#Plot time for each data set
time_breaks <- c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,300000)
crashed_y <- 800000  #height of the 'crashed' line
crashed$time_estimation <- crashed_y

results$x <- as.numeric(results$ds_name)
crashed <- crashed[order(crashed$ds_name, crashed$method),]
crashed$x <- as.numeric(crashed$ds_name) +
  unlist(lapply(split(crashed$method, crashed$ds_name),
                function(m) {
                  #skulls are moved closer together if there are many of them
                  spacing <- min(0.20, 0.55 / max(1, length(m) - 1))
                  (seq_along(m) - (length(m) + 1) / 2) * spacing
                }))

line_cols <- c("x", "method", "time_estimation")
line_crashed <- crashed[, line_cols]
line_crashed$time_estimation <- crashed_y / 1.5
line_data <- rbind(results[, line_cols], line_crashed)
line_data <- line_data[order(line_data$method, line_data$x),]

# The skull symbol has no bold face (using fontface="bold" falls back to a grey
# emoji font). It is thus drawn several times with small offsets to make it bolder.
crashed_bold <- do.call(rbind, lapply(
  list(c(0, 1), c(0.004, 1), c(-0.004, 1), c(0, 1.012), c(0, 1/1.012)),
  function(o) {
    d <- crashed
    d$x <- d$x + o[1]
    d$time_estimation <- d$time_estimation * o[2]
    d
  }))

p_time <- ggplot(results, aes(x=x, y=time_estimation, group = method)) +
  geom_hline(yintercept = 450000, linetype = "dotted", color = "grey40") +
  geom_line(data = line_data, aes(color=method), linewidth=0.5, linetype="dashed") +
  geom_point(aes(color=method, shape=method), size=3, stroke=1.2) + theme_bw(base_size = 13.5) +
  geom_text(data = crashed_bold, aes(x=x, y=time_estimation, color=method, label="☠"),
            family="Segoe UI Symbol", size=5, show.legend=FALSE) +
  scale_x_continuous(breaks = seq_along(data_sets), labels = ds_labels[data_sets],
                     expand = expansion(add = 0.55)) +
  scale_y_continuous(trans = "log1p", breaks = c(time_breaks, crashed_y),
                     labels = c(y_labels(time_breaks), "crash / NA"),
                     limits = c(NA, 950000)) +
  theme(legend.position = "top", legend.title=element_blank()) +
  xlab("Data set") + ylab("Time (s)") +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
p_time
save_plot(p_time, "real_world_applications_time")

#Plot average relative difference to the fastest model for each model
results <- results %>% group_by(ds_name) %>% mutate(min_time = min(time_estimation), 
                                                    rel_time = (time_estimation - min_time)/min_time)
mean_rel_diff <- results %>% group_by(method) %>% summarize(mean_rel_time = mean(rel_time))

p_rel_diff <- ggplot(mean_rel_diff, aes(x=method, y=mean_rel_time)) +
  geom_point(aes(color=method, shape=method), size=3, stroke=1.2) + theme_bw(base_size = 13.5) +
  theme(legend.position="none") +
  scale_y_continuous(trans = "log1p", breaks = c(0, time_breaks), labels = y_labels) +
  xlab("") + ylab("Average relative difference") +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
p_rel_diff
save_plot(p_rel_diff, "real_world_applications_rel_diff")

###Tables: runtimes, variance parameters, and coefficients######################
# LaTeX tables with the estimation results. Methods without a result are
# reported as 'crashed' (aborted with an error or ran out of memory) or as 'NA'
# (ran through, but without producing a usable fit), and methods that did not
# converge are marked with a dagger.

table_data_sets <- c("cars", "chicago_building_permits", "instEval", "MovieLens_32m",
                     "amazon_employee_access", "KDDCup09_upselling")
table_methods <- c("Cholesky (GPBoost)", "Krylov (GPBoost)", "glmmTMB", "lme4",
                   "MixedModels.jl", "R-INLA (EB)")

# Why a method has no result, as documented in the individual experiment scripts
# (see the comments next to `krylovgmms.methods` there).
failure_status <- list(
  chicago_building_permits = c("MixedModels.jl" = "crashed"),
  MovieLens_32m = c("Cholesky (GPBoost)" = "crashed", "glmmTMB" = "crashed",
                    "MixedModels.jl" = "crashed", "R-INLA (EB)" = "crashed"),
  KDDCup09_upselling = c("glmmTMB" = "NA", "MixedModels.jl" = "NA",
                         "R-INLA (EB)" = "crashed")
)

failure_label <- function(d, method) {
  status <- failure_status[[d]][method]
  if (is.na(status)) {
    warning("No documented failure reason for ", method, " on ", d,
            "; reporting it as 'NA'.", call. = FALSE)
    return("NA")
  }
  unname(status)
}

# A method that stopped at a clearly worse optimum than the best one found for
# the same data set did not converge properly. R-INLA is excluded because it
# integrates the fixed effects out and its criterion is therefore not on the
# same scale as the ML values of the other methods (see experiment_helpers.R).
nll_tolerance <- 1
comparable_nll_methods <- setdiff(table_methods, "R-INLA (EB)")

# Appended to every caption to explain the entries for methods without a result
# and the marker for methods that did not converge.
table_note <- paste0(
  "'crashed' means that estimation aborted with an error or ran out of memory, ",
  "and 'NA' that estimation terminated but returned missing values for the ",
  "parameter estimates or the negative log-marginal likelihood at the optimum. ",
  "A '*' in front of a method marks fits whose negative log-likelihood at the ",
  "optimum is more than ", nll_tolerance, " above the smallest one obtained for ",
  "the data set, i.e. that did not converge properly."
)

# Explains the dagger on the objective values that are not comparable to the ML
# ones. Only the runtime table reports these values; the other two captions just
# state why the convergence criterion does not apply to all methods.
nll_note <- paste0(
  "This convergence criterion applies only to the maximum ",
  "marginal-likelihood-based implementations. The values marked with ",
  "'$\\dagger$' are the negative log-marginal likelihoods reported by ",
  "\\texttt{R-INLA}. They are obtained with the fixed effects integrated out ",
  "and are thus restricted (REML-type) criteria evaluated with the INLA ",
  "approximation, see Section \\ref{exp_setting}. They are not on the same ",
  "scale as the maximum likelihood values of the other implementations and are ",
  "therefore not comparable to them."
)
estimates_note <- paste0(
  "This convergence criterion applies only to the maximum ",
  "marginal-likelihood-based implementations since the objective function ",
  "reported by \\texttt{R-INLA} is not on the same scale, see Table ",
  "\\ref{table:real_world_times}."
)
not_converged <- function(r) {
  comparable <- r$method %in% comparable_nll_methods & is.finite(r$nll_optimum)
  if (!any(comparable)) return(rep(FALSE, nrow(r)))
  comparable & r$nll_optimum > min(r$nll_optimum[comparable]) + nll_tolerance
}

latex_escape <- function(x) gsub("_", "\\\\_", x)

#value in the form $a.bc \times 10^{d}$
sci_label <- function(x) {
  if (length(x) != 1 || is.na(x)) return("")
  if (x == 0) return("$0.00 \\times 10^{0}$")
  e <- floor(log10(abs(x)))
  m <- x / 10^e
  if (abs(round(m, 2)) >= 10) {  #e.g. 9.999 -> 1.00 x 10^(e+1)
    m <- m / 10
    e <- e + 1
  }
  sprintf("$%.2f \\times 10^{%d}$", m, e)
}

num_label <- function(x, digits = 1) {
  if (length(x) != 1 || is.na(x)) return("")
  formatC(x, format = "f", digits = digits)
}

# The objective of a method that is not on the same scale as the ML values of
# the other methods (R-INLA, see above) is marked with a dagger.
nll_cell <- function(row) {
  value <- num_label(row$nll_optimum)
  if (!nzchar(value) || row$method %in% comparable_nll_methods) return(value)
  paste0("$", value, "^{\\dagger}$")
}

#value of column 'name', NA if the column does not exist for this data set
col_value <- function(row, name) if (name %in% names(row)) row[[name]] else NA

get_table_results <- function(d) {
  r <- read_results(d)$all_results
  r[match(table_methods, r$method),]
}

# 'columns' is a list of functions that map a row of the results to a table cell
make_latex_table <- function(columns, header, align, caption, label) {
  lines <- c("\\begin{table}[ht!]", "\\centering",
             paste0("\\begin{tabular}{", align, "}"), "  \\hline",
             paste0("data set & method & ", header, " \\\\"), "  \\hline")
  for (d in table_data_sets) {
    r <- get_table_results(d)
    diverged <- not_converged(r)
    for (i in seq_along(table_methods)) {
      if (is.na(r$time_estimation[i])) {  #crashed or no usable fit
        cells <- c(failure_label(d, table_methods[i]), rep("", length(columns) - 1))
      } else {
        cells <- sapply(columns, function(f) f(r[i,]))
      }
      method_label <- paste0(if (diverged[i]) "*", table_methods[i])
      lines <- c(lines, paste0("  ", if (i == 1) latex_escape(ds_labels[[d]]) else "",
                               " & ", method_label, " & ",
                               paste(cells, collapse = " & "), " \\\\"))
    }
    lines <- c(lines, "  \\hline")
  }
  cat(c(lines, "\\end{tabular}",
        paste0("\\caption{", caption, "}"),
        paste0("\\label{", label, "}"), "\\end{table}", ""), sep = "\n")
}

#Table: runtimes and negative log-likelihood at the optimum
make_latex_table(
  columns = list(function(x) num_label(x$time_estimation), nll_cell),
  header = "Time (s) & nll\\_optimum", align = "llrr",
  caption = paste("Time for parameter estimation and negative log-marginal likelihood",
                  "at the optimum for different real-world data sets and models.",
                  table_note, nll_note),
  label = "table:real_world_times")

#Table: variance parameters
make_latex_table(
  columns = lapply(c("sigma", "sigma2_1", "sigma2_2", "sigma2_3"),
                   function(nm) { force(nm); function(x) sci_label(col_value(x, nm)) }),
  header = "$\\sigma^2$ & $\\sigma_1^2$ & $\\sigma_2^2$ & $\\sigma_3^2$", align = "llcccc",
  caption = paste("Estimates for $\\sigma^2$, $\\sigma_1^2$, $\\sigma_2^2$, and $\\sigma_3^2$",
                  "for different real-world data sets and models. For reasons of space,",
                  "we do not report estimates for other variance parameters.",
                  table_note, estimates_note),
  label = "table:real_world_variances")

#Table: regression coefficients
make_latex_table(
  columns = lapply(c("beta_0", "beta_1", "beta_2", "beta_3"),
                   function(nm) { force(nm); function(x) sci_label(col_value(x, nm)) }),
  header = "$\\beta_0$ & $\\beta_1$ & $\\beta_2$ & $\\beta_3$", align = "llcccc",
  caption = paste("Estimates for $\\beta_0$, $\\beta_1$, $\\beta_2$, and $\\beta_3$",
                  "for different real-world data sets and models. For reasons of space,",
                  "we do not report estimates for other coefficients.",
                  table_note, estimates_note),
  label = "table:real_world_coefficients")
