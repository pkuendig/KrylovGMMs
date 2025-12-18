################################################################################
# Preconditioner comparison 
################################################################################
library(gpboost)
library(ggplot2)
library(grid)
library(dplyr)
library(patchwork)

## Load function for simulating data and set some settings
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")
set.seed(1)

PRECONDITIONER_plot_names <- list()
PRECONDITIONER_plot_names[["symmetric_successive_over_relaxation"]] <- "SSOR"
PRECONDITIONER_plot_names[["zero_infill_incomplete_cholesky"]] <- "ZIC"
PRECONDITIONER_plot_names[["diagonal"]] <- "Diagonal"
PRECONDITIONER_plot_names[["none"]] <- "None"

run_simulations <- TRUE # If TRUE, calculations are run, otherwise only results are loaded and plots are generated
run_cholesky <- TRUE # If TRUE, calculations for Cholesky are run (slow), otherwise only results are loaded
n_rep <- 100 # number of simulation runs

## Function for creating plots
make_plots <- function(Itresults_1, CHOLtime=NULL, CHOL_ll=NULL, 
                       numeric_x_scale = FALSE, chol_results = NULL, title = NULL,
                       adj_size = 1.5) {
  
  if (numeric_x_scale) {
    Itresults_1[[2]] <- as.numeric(as.character(Itresults_1[[2]]))
  }
  xvar <- names(Itresults_1)[2]
  xvar_sym <- rlang::sym(xvar)
  
  p1 <- ggplot(Itresults_1, aes(x=!!xvar_sym, y=SD, color=preconditioner, shape=preconditioner, group=preconditioner)) +
    geom_point(size = 2) + geom_line(linewidth = 1) + 
    scale_y_continuous(trans = "log1p", breaks = c(0.3,1,2,5,10,20,50,100,200,500,1000)) +
    scale_color_brewer(type = "qual", palette = 6, name = "Preconditioner") + scale_shape_discrete(name = "Preconditioner") + 
    labs(color = "") + ylab("Std. dev. of log-likelihood") + xlab("") + theme_bw() + 
    theme(legend.position = "top", legend.text = element_text(size = 12*adj_size), legend.title = element_text(size = 12*adj_size),
          plot.margin = margin(t = 0, r = 40, b = 0, l = 5),
          axis.title.y = element_text(size = 12*adj_size),
          axis.text.x = element_text(size = 12*adj_size), axis.text.y = element_text(size = 12*adj_size)) + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, margin = margin(t = 5), size = 16*adj_size, face = "bold"))
  }
  
  if (max(Itresults_1$CG_iter, na.rm = TRUE) <= 10){
    cg_breaks <- seq(0,10,2)
  } else if (max(Itresults_1$CG_iter, na.rm = TRUE) <= 20){
    cg_breaks <- seq(0,20,4)
  } else if (max(Itresults_1$CG_iter, na.rm = TRUE) <= 30){
    cg_breaks <- seq(0,30,5)
  } else if (max(Itresults_1$CG_iter, na.rm = TRUE) <= 100){
    cg_breaks <- seq(10,100,10)
  } else {
    cg_breaks <- seq(50,1000,50)
  }
  cg_breaks <- cg_breaks[cg_breaks <= max(Itresults_1$CG_iter, na.rm = TRUE)]
  scale_factor <- max(Itresults_1$time, na.rm = TRUE) / max(Itresults_1$CG_iter, na.rm = TRUE)
  cg_lines  <- cg_breaks * scale_factor  # transform to time scale
  p2 <- ggplot(Itresults_1, aes(x=!!xvar_sym, y=time, color=preconditioner, shape=preconditioner, group=preconditioner)) +
    geom_point(aes(y = time), size = 2) + geom_line(aes(y = time, linetype = "Time (left axis)"), linewidth = 1) +
    scale_color_brewer(type = "qual", palette=6) + labs(color = "") + ylab("Time (s)") + 
    theme_bw() + theme(axis.title.y.right = element_text(margin = margin(l = 10, r = -60), size = 12*adj_size), 
                       axis.title.y = element_text(size = 12*adj_size), axis.title.x = element_text(size = 14*adj_size),
                       axis.text.x = element_text(size = 12*adj_size), axis.text.y = element_text(size = 12*adj_size), 
                       plot.margin = margin(t = 0, r = 40, b = 0, l = 5),
                       legend.position = "bottom", legend.text = element_text(size = 11*adj_size), legend.key.width = unit(2, "cm"),
                       legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0))  + 
    guides(color = "none", shape = "none") + 
    geom_point(aes(y = CG_iter * scale_factor), size = 2) +
    geom_line(aes(y = CG_iter * scale_factor, linetype = "CG iterations (right axis)"), linewidth = 1) +
    geom_hline(yintercept = cg_lines, color = "gray80", linetype = "dotted") +
    scale_y_continuous( name = "Time (s)", sec.axis = sec_axis(trans = ~ . / scale_factor,
                                                               name = "CG iterations", breaks = cg_breaks)) +
    scale_linetype_manual(name = "", values = c("Time (left axis)" = "solid",
                                                "CG iterations (right axis)" = "dashed"), limits = c("Time (left axis)", "CG iterations (right axis)"))
  
  if (!is.null(chol_results)) {
    p2 <- p2 + geom_point(data = chol_results, aes(x = !!xvar_sym, y = time, group = 1), 
                          inherit.aes = FALSE, size = 3*adj_size, shape = 17 ) + 
      annotate("text", x = Inf, y = Inf, label = "▲  Cholesky (time)", hjust = 1.1, vjust = 2, size = 4)
  }
  
  if (numeric_x_scale) {
    p1 <- p1 + scale_x_continuous(breaks = unique(Itresults_1[[2]]))
    p2 <- p2 + scale_x_continuous(breaks = unique(Itresults_1[[2]]))
  }
  if (!is.null(CHOLtime)) {
    p2 <- p2 + annotate("text", label = paste0("Cholesky: ", round(CHOLtime, digits = 1), " (s)"), 
                        x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size=5*adj_size)
  }
  if (!is.null(CHOL_ll)) {
    p1 <- p1 + annotate("text", label = paste0("Negative log-lik. (Cholesky): ", round(CHOL_ll, digits = 1)), 
                        x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size=5*adj_size)
  }
  combined <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first")
  return(combined)
}

################################################################################
# Varying the number of probe vectors for the SLQ method
################################################################################

plot_height <- 10
plot_width <- 8

m1 <- 50000
n <- (2*m1) * 10

NUM_RAND_VEC_TRACE <- c(10, 20, 50, 100)
PRECONDITIONER <- c("symmetric_successive_over_relaxation", "zero_infill_incomplete_cholesky", "diagonal", "none")

for (snr in c(2,10)) {
  sigma2_1 <- 0.5^2
  sigma2_2 <- 0.5^2
  sigma2 <- (sigma2_2 + sigma2_1) / snr
  
  for (likelihood in c("gaussian", "bernoulli_logit")) {
    if (snr == 2 | likelihood == "gaussian") {
      if (likelihood == "gaussian"){
        cov_pars <- c(sigma2,sigma2_1,sigma2_2)
      } else {
        cov_pars <- c(sigma2_1,sigma2_2)
      }
      for (balanced in c(TRUE, FALSE)) {
        cat("\n\n")
        print(paste0("******** Starting ll vs. t, likelihood = ", likelihood, 
                     ", snr = ", snr, ", balanced = ",balanced))
        
        file_name <- paste0("compare_P_nll_",likelihood,"_vs_t")
        if (snr != 2) file_name <- paste0(file_name,"_high_SNR")
        if(balanced) file_name <- paste0(file_name,"_balanced")
        set.seed(m1)
        the_data <- make_data(n=n, m1=m1, sigma2=sigma2, sigma2_1=sigma2_1, sigma2_2=sigma2_2,
                              likelihood=likelihood, balanced=balanced)
        if (run_simulations) {
          Itresults <- data.frame(matrix(nrow=length(NUM_RAND_VEC_TRACE)*length(PRECONDITIONER)*n_rep,ncol = 5))
          colnames(Itresults) <- c("preconditioner", "t", "negLL", "time", "CG_iter")
          i = 1
          for(t in 1:length(NUM_RAND_VEC_TRACE)){
            cat("\n")
            print(paste0("num_rand_vec_trace = ", NUM_RAND_VEC_TRACE[t]))
            for(p in 1:length(PRECONDITIONER)){
              print(paste0("cg_preconditioner_type = ", PRECONDITIONER[p]))
              for(r in 1:n_rep){
                cat(".")
                if(r%%10 == 0) cat(paste0("r=", r))
                Itmodel <- GPModel(group_data = the_data$group_data,
                                   likelihood=likelihood, matrix_inversion_method = "iterative")
                Itmodel$set_optim_params(params = list(init_cov_pars=cov_pars,
                                                       num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                                       cg_preconditioner_type = PRECONDITIONER[p],
                                                       seed_rand_vec_trace=i))
                Itresults$preconditioner[i] <- PRECONDITIONER_plot_names[[PRECONDITIONER[p]]]
                Itresults$t[i] <- NUM_RAND_VEC_TRACE[t]
                Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=cov_pars, y=the_data$y))[3]
                Itresults$CG_iter[i] <- Itmodel$get_num_cg_steps()
                i = i+1
                gc()
              }
              cat("\n")
            }
          }
          Itresults$preconditioner <- factor(Itresults$preconditioner, levels=unlist(PRECONDITIONER_plot_names, use.names = FALSE))
          Itresults$t <- as.factor(Itresults$t)
          saveRDS(Itresults, file=paste0("./../results/", file_name,".rds"))
        } else {
          Itresults <- readRDS(file=paste0("./../results/", file_name,".rds"))
        }
        
        ###Cholesky#####################################################################
        if (snr == 2) {
          if(run_cholesky){
            print("Starting cholesky-based calculations ")
            cat("\n")
            CHOLmodel <- GPModel(group_data = the_data$group_data, likelihood=likelihood, matrix_inversion_method = "cholesky")
            CHOLmodel$set_optim_params(params = list(init_cov_pars=cov_pars))
            CHOLtime <- system.time(CHOLresult <- CHOLmodel$neg_log_likelihood(cov_pars=cov_pars, y=the_data$y))[3]
            chol_results <- c(CHOLresult, CHOLtime)
            saveRDS(chol_results, file=paste0("./../results/", file_name,"_chol.rds"))
          } else {
            chol_results <- readRDS(file=paste0("./../results/", file_name,"_chol.rds"))
            CHOLresult <- chol_results[1]
            CHOLtime <- chol_results[2]
          }
        } else {
          CHOLresult <- NULL
          CHOLtime <- NULL
        }
        ## Plot results
        if (balanced) {
          title <- "Balanced"
        } else {
          title <- "Unbalanced"
        }
        if (snr == 2) { 
          adj_size <- 1.5
          p1 <- ggplot(Itresults, aes(x=t, y=negLL, fill=preconditioner)) + 
            geom_hline(yintercept=CHOLresult, linetype = "dashed") + theme_bw() +
            geom_boxplot() + labs(fill  = "") + ylab("log-likelihood") +
            scale_fill_brewer(type = "qual", palette=6, labels = scales::parse_format()) + 
            scale_y_continuous(n.breaks=8) +
            guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + 
            theme(axis.title.x=element_blank(),  axis.text.x=element_blank(), 
                  axis.ticks.x=element_blank(),  legend.position = "top", 
                  axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0), size = 12*adj_size),
                  axis.text.y = element_text(size = 12*adj_size),
                  legend.text = element_text(size = 12*adj_size), legend.title = element_text(size = 12*adj_size))+ 
            ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, margin = margin(t = 5), size = 16*adj_size, face = "bold"))
          p2 <- ggplot(Itresults, aes(x=t, y=time, color=preconditioner, shape=preconditioner)) +
            stat_summary(aes(group = preconditioner), fun = mean, geom = 'line', linewidth=1, alpha=0.9) + 
            stat_summary(aes(group = preconditioner), fun = mean, geom = 'point', size=2) +
            scale_color_brewer(type = "qual", palette=6) + labs(color = "") + ylab("Time (s)") +
            theme_bw() + theme(legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0), size = 12*adj_size),
                               axis.title.x = element_text(size = 14*adj_size), legend.text = element_text(size = 11*adj_size),
                               axis.text.x = element_text(size = 12*adj_size), axis.text.y = element_text(size = 12*adj_size)) +
            annotate("text", label = paste0("Cholesky: ", round(CHOLtime, digits = 1), "s"), x = 4, y = 0.3, size=4*adj_size) 
          combined <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first")
          grid.newpage()
          grid.draw(combined)
          ggsave(combined, file=paste0("./../plots/", file_name,"_hist.png"),height=plot_height,width=plot_width)
        }
        Itresults_1 <- Itresults %>% group_by(preconditioner,t) %>% summarise(SD = sd(negLL), time = mean(time), CG_iter = mean(CG_iter), .groups = "drop")
        combined <- make_plots(Itresults_1, CHOLtime=CHOLtime, CHOL_ll = CHOLresult, numeric_x_scale=TRUE, title=title, adj_size=1.5)
        grid.newpage()
        grid.draw(combined)
        ggsave(combined, file=paste0("./../plots/", file_name,".png"),height=plot_height,width=plot_width)
      
      }
    }
  }
}

################################################################################
# Varying the random effects dimension and the sample size
################################################################################

sigma2 <- 0.5^2
sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
cov_pars <- c(sigma2,sigma2_1,sigma2_2)

M1s <- c(25000,50000,100000,250000)
PRECONDITIONER <- c("symmetric_successive_over_relaxation", "zero_infill_incomplete_cholesky", "diagonal", "none")

cat("\n\n")
print("******** Starting ll vs. m ")
file_name <- "compare_P_nll_Gaussian_vs_m"
if (run_simulations) {
  Itresults <- data.frame(matrix(nrow=length(M1s)*length(PRECONDITIONER)*n_rep,ncol = 5))
  colnames(Itresults) <- c("preconditioner", "m", "negLL", "time", "CG_iter")
  i = 1
  for(m1 in M1s){
    n <- (2*m1) * 10
    set.seed(m1)
    the_data <- make_data(n=n, m1=m1, sigma2=sigma2, sigma2_1=sigma2_1, sigma2_2=sigma2_2,
                          likelihood="gaussian", balanced = FALSE)
    cat("\n")
    print(paste0("m1 = ", m1))
    for(p in 1:length(PRECONDITIONER)){
      print(paste0("cg_preconditioner_type = ", PRECONDITIONER[p]))
      for(r in 1:n_rep){
        cat(".")
        if(r%%10 == 0) cat(paste0("r=", r))
        Itmodel <- GPModel(group_data = the_data$group_data,
                           likelihood="gaussian", matrix_inversion_method = "iterative")
        Itmodel$set_optim_params(params = list(init_cov_pars=cov_pars,
                                               num_rand_vec_trace=50,
                                               cg_preconditioner_type = PRECONDITIONER[p],
                                               seed_rand_vec_trace=i))
        Itresults$preconditioner[i] <- PRECONDITIONER_plot_names[[PRECONDITIONER[p]]]
        Itresults$m[i] <- 2*m1
        Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=cov_pars, y=the_data$y))[3]
        Itresults$CG_iter[i] <- Itmodel$get_num_cg_steps()
        i = i+1
        gc()
      }
      cat("\n")
    }
  }
  Itresults$preconditioner <- factor(Itresults$preconditioner, levels=unlist(PRECONDITIONER_plot_names, use.names = FALSE))
  Itresults$m <- as.factor(Itresults$m)
  saveRDS(Itresults, file=paste0("./../results/", file_name,".rds"))
} else {
  Itresults <- readRDS(file=paste0("./../results/", file_name,".rds"))
}
## Plot results
Itresults_1 <- Itresults %>% group_by(preconditioner,m) %>% summarise(SD = sd(negLL), time = mean(time), CG_iter = mean(CG_iter), .groups = "drop")
combined <- make_plots(Itresults_1, title ="Varying m")
grid.newpage()
grid.draw(combined)
ggsave(combined, file=paste0("./../plots/", file_name,".png"),height=plot_height,width=plot_width)


################################################################################
# Varying the regularity of the random effects design
################################################################################

sigma2 <- 0.5^2
sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
cov_pars <- c(sigma2,sigma2_1,sigma2_2)

OverDisps <- c(0.1,1,10,100)
PRECONDITIONER <- c("symmetric_successive_over_relaxation", "zero_infill_incomplete_cholesky", "diagonal", "none")
m1 <- 50000
n <- (2*m1) * 10

cat("\n\n")
print("******** Starting ll vs. over_disp")
file_name <- "compare_P_nll_Gaussian_vs_over_disp"
if (run_simulations) {
  Itresults <- data.frame(matrix(nrow=length(OverDisps)*length(PRECONDITIONER)*n_rep,ncol = 5))
  colnames(Itresults) <- c("preconditioner", "over_disp", "negLL", "time", "CG_iter")
  i = 1
  for(over_disp in OverDisps){
    set.seed(over_disp)
    the_data <- make_data(n=n, m1=m1, sigma2=sigma2, sigma2_1=sigma2_1, sigma2_2=sigma2_2,
                          likelihood="gaussian", size_neg_bin = 1/over_disp, balanced = FALSE)
    cat("\n")
    print(paste0("over_disp = ", over_disp))
    for(p in 1:length(PRECONDITIONER)){
      print(paste0("cg_preconditioner_type = ", PRECONDITIONER[p]))
      for(r in 1:n_rep){
        cat(".")
        if(r%%10 == 0) cat(paste0("r=", r))
        Itmodel <- GPModel(group_data = the_data$group_data,
                           likelihood="gaussian", matrix_inversion_method = "iterative")
        Itmodel$set_optim_params(params = list(init_cov_pars=cov_pars,
                                               num_rand_vec_trace=50,
                                               cg_preconditioner_type = PRECONDITIONER[p],
                                               seed_rand_vec_trace=i))
        Itresults$preconditioner[i] <- PRECONDITIONER_plot_names[[PRECONDITIONER[p]]]
        Itresults$over_disp[i] <- over_disp
        Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=cov_pars, y=the_data$y))[3]
        Itresults$CG_iter[i] <- Itmodel$get_num_cg_steps()
        i = i+1
        gc()
      }
      cat("\n")
    }
  }
  Itresults$preconditioner <- factor(Itresults$preconditioner, levels=unlist(PRECONDITIONER_plot_names, use.names = FALSE))
  Itresults$over_disp <- as.factor(Itresults$over_disp)
  saveRDS(Itresults, file=paste0("./../results/", file_name,".rds"))
} else {
  Itresults <- readRDS(file=paste0("./../results/", file_name,".rds"))
}
## Plot results
Itresults_1 <- Itresults %>% group_by(preconditioner,over_disp) %>% summarise(SD = sd(negLL), time = mean(time), CG_iter = mean(CG_iter), .groups = "drop")
combined <- make_plots(Itresults_1, title = "Varying irregularity of the design")
grid.newpage()
grid.draw(combined)
ggsave(combined, file=paste0("./../plots/", file_name,".png"),height=plot_height,width=plot_width)

