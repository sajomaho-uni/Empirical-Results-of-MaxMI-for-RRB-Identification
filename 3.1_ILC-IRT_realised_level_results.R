# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for reproducing the results of employing the MaxMI-22
# method in the ILC-IRT framework at the realised level

# Author:
# Santeri Holopainen
# Turku Research Institute for Learning Analytics (TRILA)
# University of Turku

# Date: 11.9.2026

# R version: 4.5.0

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# 1 Setup ----

# Clear the environment
rm(list = ls())

# Load the following packages
library(tidyverse)
library(tidylog)
library(crayon)
library(ggpubr)
library(ggh4x)
library(cubature)
library(MASS)
library(mixtools)

# Settings for tidylog
crayon = function(x) cat(green(x), sep = "\n")
options("tidylog.display" = list(crayon))
rm(crayon)

# 2 Form a data frame of the simulation parameters ----

# See Table 1 in the Supplement

# Sample size:
N <- c(500, 2500)

# Person parameters:

# Mean and variance of the latent engagement tendency (phi)
mu_phi <- 0; var_phi <- 3.5

# Mean and variance of the latent ability (theta)
mu_theta <- 0; var_theta <- 1

# Mean and variance of the latent speed (tau)
mu_tau <- 0; var_tau <- 0.05

# Correlation between phi and theta
cor_phi_theta <- .55

# Correlation between phi and tau
cor_phi_tau <- .2

# Correlations between theta and tau
cor_theta_tau <- c(-0.4, 0.4)

# Item parameters:

# Item engagement difficulty
iota <- c(-4.27, -3.29)

# Chance level
g <- 1/4

# Item difficulty
b <- c(-2, -1, 0, 1, 2)

# In the study, we use the Rasch model for the probability of answering correctly in the SB component
# So, set item discrimination to 1 and pseudo-guessing to 0:
a <- 1
c <- 0

# Mean of the disengaged log RTs
mu_c <- 3

# SD of the disengaged log RTs
sigma_c <- sqrt(1.95)

# Time intensity
beta <- c(3.25, 4.25)

# Time discrimination
alpha <- sqrt(0.15)

# Add all parameter combinations to a data frame
params_real_lvl <- expand.grid(N,
                               mu_phi, mu_theta, mu_tau, 
                               var_phi, var_theta, var_tau, 
                               cor_phi_theta, cor_phi_tau, cor_theta_tau,
                               iota,
                               mu_c, sigma_c, alpha, beta,
                               g, a, b, c)

# Rename columns
colnames(params_real_lvl) <- c("N",
                               "mu_phi", "mu_theta", "mu_tau",
                               "var_phi", "var_theta", "var_tau",
                               "cor_phi_theta", "cor_phi_tau", "cor_theta_tau",
                               "iota",
                               "mu_c", "sigma_c", "alpha", "beta",
                               "g", "a", "b", "c")

# Compute also the covariances between the person parameters
params_real_lvl <-
  params_real_lvl %>%
  mutate(cov_phi_theta = cor_phi_theta*sqrt(var_phi*var_theta),
         cov_phi_tau = cor_phi_tau*sqrt(var_phi*var_tau),
         cov_theta_tau = cor_theta_tau*sqrt(var_theta*var_tau))

# Arrange the data frame and add simulation condition number
params_real_lvl <-
  params_real_lvl %>%
  arrange(N, cor_theta_tau, iota, beta, b) %>%
  mutate(sim = 1:n()) %>%
  relocate(sim, .before = N)

# save
save(params_real_lvl, file = "params_real_lvl.rda")

# 3 Compute the conditional marginal probabilities of a correct response for SB test-takers at the population level ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_real_lvl.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "p_u1_cond_d1_f" over the 80 different conditions
res <- lapply(params_real_lvl$sim, function(x) p_u1_cond_d1_f(k = x, params = params_real_lvl))

# When finished, combine to a data frame
p_u1_cond_d1_real_lvl_df <- bind_rows(lapply(res, function(x) data.frame(p_u1_cond_d1 = x)), .id = "sim")

# save
save(p_u1_cond_d1_real_lvl_df, file = "p_u1_cond_d1_real_lvl_df.rda")

# 4 Simulation ----

# . 4.1 Generate the seed numbers for data generation and threshold computation ----

# Clear the environment
rm(list = ls())

# Generate seed numbers
seeds <- expand.grid(1:100, 1:80)
colnames(seeds) <- c("rep", "sim")
seeds <- 
  seeds %>% 
  mutate(seed_data = 1:n(), # seeds for data generation
         seed_thr = (n()+1):(2*n())) # seeds for threshold computation

# Save the seed data frame
save(seeds, file = "seeds_for_simulations.rda")

# . 4.2 Generate data ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_real_lvl.rda")

# Import the seeds
load("seeds_for_simulations.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "gen_data" over all repetitions and simulation conditions
system.time(
    sim_data <- lapply(1:100, # repetitions
                        function(x){
                          lapply(1:80, # simulation conditions
                                 function(y){
                                   # Print repetition number and simulation condition number for keeping track of the progress
                                   print(paste0("REP: ", x, ", SIM: ", y))
                                 
                                   # Get seed number
                                   seed_k <- seeds$seed_data[seeds$rep == x & seeds$sim == y]
                                 
                                   # Generate the data
                                   gen_data(k = y, params = params_real_lvl, seed = seed_k)
                               }
                        )
                      }
  )
)

# Save the generated data
save(sim_data, file = "sim_data.rda")

# . 4.3 Compute posterior probabilities of belonging to the engaged component P(Delta = 1 | l, y) ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Pre-compute P(Delta = 1 | iota) for all iota
integrand <- function(phi, iota, mu_phi, sigma_phi) p_engagement(phi = phi, iota = iota)*dnorm(phi, mean = mu_phi, sd = sigma_phi)
p_d1_df <- data.frame(iota = unique(params_real_lvl$iota), p_d1 = NA)
p_d1_df$p_d1[1] <- integrate(integrand, 
                             lower = -100, upper = 100,
                             iota = p_d1_df$iota[1],
                             mu_phi = 0, 
                             sigma_phi = sqrt(3.5))$value
p_d1_df$p_d1[2] <- integrate(integrand, 
                             lower = -100, upper = 100,
                             iota = p_d1_df$iota[2],
                             mu_phi = 0, 
                             sigma_phi = sqrt(3.5))$value

# Apply "compute_p_d1_posterior" over all repetitions and simulation conditions
system.time(
  p_d1_cond_ly_real_lvl <- lapply(1:100, # repetitions
                                     function(x){
                                       lapply(1:80, # simulation conditions
                                              function(y){
                                                # Print repetition number and simulation condition number for keeping track of the progress
                                                print(paste0("REP: ", x, ", SIM: ", y))
                                               
                                                # Get data
                                                temp_data <- sim_data[[x]][[y]]
                                               
                                                # Compute MR posterior
                                                compute_p_d1_posterior(k = y, 
                                                                       params = params_real_lvl, 
                                                                       sim_data = temp_data, 
                                                                       n_gh = 25, 
                                                                       p_d1_df = p_d1_df)
                                                }
                                              )
                                       }
                                     )
  )

# Save the results
save(p_d1_cond_ly_real_lvl, file = "p_d1_cond_ly_real_lvl.rda")

# . 4.4 Compute thresholds ----

# . . 4.4.1 Compute population-level MaxMI-22 thresholds (needed in assessing bias of MaxMI-22) ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_real_lvl.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_pop_lvl_MaxMI22_f" over the 80 different conditions
thresholds_pop_lvl_MaxMI22_for_simulations <- lapply(params_real_lvl$sim, function(x) thr_pop_lvl_MaxMI22_f(sim = x, params = params_real_lvl))

# Save the thresholds
save(thresholds_pop_lvl_MaxMI22_for_simulations, file = "thresholds_pop_lvl_MaxMI22_for_simulations.rda")

# . . 4.4.2 Compute thresholds by minimising the misclassification rate based on the posterior probability P(Delta = 1 | l, y) ----

# Clear the environment
rm(list = ls())

# Import the posterior probs
load("p_d1_cond_ly_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_real_lvl_MinMR_f" over all repetitions and simulation conditions
thresholds_real_lvl_MinMR <- 
  lapply(1:100, # repetitions
         function(x){
           lapply(1:80, # simulation conditions
                  function(y){
                    # Print repetition number and simulation condition number for keeping track of the progress
                    print(paste0("REP: ", x, ", SIM: ", y))
                    
                    # Get data
                    temp_data <- sim_data[[x]][[y]]
                    
                    # Get the posterior probabilities
                    temp_p_d1_cond_ly <- p_d1_cond_ly_real_lvl[[x]][[y]]
                    
                    # Find the threshold(s)
                    thr_real_lvl_MinMR_f(sim_data_k = temp_data, p_d1_cond_ly_k = temp_p_d1_cond_ly)
                  })
         })

# Save
save(thresholds_real_lvl_MinMR, file = "thresholds_real_lvl_MinMR.rda")

# . . 4.4.3 Compute MaxMI-22 thresholds ----

# Clear the environment
rm(list = ls())

# Import the posterior probs
load("p_d1_cond_ly_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_real_lvl_MaxMI22_f" over all repetitions and simulation conditions
thresholds_real_lvl_MaxMI22 <- 
  lapply(1:100, # repetitions
         function(x){
           lapply(1:80, # simulation conditions
                  function(y){
                    # Print repetition number and simulation condition number for keeping track of the progress
                    print(paste0("REP: ", x, ", SIM: ", y))
                    
                    # Get data
                    temp_data <- sim_data[[x]][[y]]
                    
                    # Get the posterior probabilities
                    temp_p_d1_cond_ly <- p_d1_cond_ly_real_lvl[[x]][[y]]
                    
                    # Find the threshold(s)
                    thr_real_lvl_MaxMI22_f(sim_data_k = temp_data, 
                                           p_d1_cond_ly_k = temp_p_d1_cond_ly,
                                           n_cut = 10, 
                                           limit = T)
                  })
         })

# Save
save(thresholds_real_lvl_MaxMI22, file = "thresholds_real_lvl_MaxMI22.rda")

# . . 4.4.4 Compute MaxMI-22 thresholds with smoothing ----

# Clear the environment
rm(list = ls())

# Import the posterior probs
load("p_d1_cond_ly_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_real_lvl_MaxMI22_smooth_f" over all repetitions and simulation conditions
thresholds_real_lvl_MaxMI22_smooth <- 
  lapply(1:100, # repetitions
         function(x){
           lapply(1:80, # simulation conditions
                  function(y){
                    # Print repetition number and simulation condition number for keeping track of the progress
                    print(paste0("REP: ", x, ", SIM: ", y))
                    
                    # Get data
                    temp_data <- sim_data[[x]][[y]]
                    
                    # Get the posterior probabilities
                    temp_p_d1_cond_ly <- p_d1_cond_ly_real_lvl[[x]][[y]]
                    
                    # Find the threshold(s)
                    thr_real_lvl_MaxMI22_smooth_f(sim_data_k = temp_data, 
                                                  p_d1_cond_ly_k = temp_p_d1_cond_ly,
                                                  n_cut = 10, 
                                                  limit = T)
                  })
         })

# Save
save(thresholds_real_lvl_MaxMI22_smooth, file = "thresholds_real_lvl_MaxMI22_smooth.rda")

# . . 4.4.5 Compute CUMP thresholds ----

# Clear the environment
rm(list = ls())

# Import the posterior probs
load("p_d1_cond_ly_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_real_lvl_CUMP_f" over all repetitions and simulation conditions
thresholds_real_lvl_CUMP <- 
  lapply(1:100, # repetitions
         function(x){
           lapply(1:80, # simulation conditions
                  function(y){
                    # Print repetition number and simulation condition number for keeping track of the progress
                    print(paste0("REP: ", x, ", SIM: ", y))
                    
                    # Get data
                    temp_data <- sim_data[[x]][[y]]
                    
                    # Get the posterior probabilities
                    temp_p_d1_cond_ly <- p_d1_cond_ly_real_lvl[[x]][[y]]
                    
                    # Find the threshold(s)
                    thr_real_lvl_CUMP_f(sim_data_k = temp_data, 
                                        p_d1_cond_ly_k = temp_p_d1_cond_ly,
                                        n_cut = 10)
                  })
         })

# Save
save(thresholds_real_lvl_CUMP, file = "thresholds_real_lvl_CUMP.rda")

# . . 4.4.6 Compute MLN thresholds ----

# Clear the environment
rm(list = ls())

# Import the posterior probs
load("p_d1_cond_ly_real_lvl.rda")

# Import the generated data
load("sim_data.rda")

# Import the seed data frame
load("seeds_for_simulations.rda")

# Read the user-defined functions from the "3.0_ILC-IRT_realised_level_functions.R" script
source("3.0_ILC-IRT_realised_level_functions.R")

# Apply "thr_real_lvl_MLN_f" over all repetitions and simulation conditions
thresholds_real_lvl_MLN <- 
  lapply(1:100, # repetitions
         function(x){ 
           lapply(1:80, # simulation conditions
                  function(y){
                    # Print repetition number and simulation condition number for keeping track of the progress
                    print(paste0("REP: ", x, ", SIM: ", y))
                    
                    # Get data
                    temp_data <- sim_data[[x]][[y]]
                    
                    # Get the posterior probabilities
                    temp_p_d1_cond_ly <- p_d1_cond_ly_real_lvl[[x]][[y]]
                    
                    # Get the seed number (for fitting a mixture model)
                    temp_seed <- seeds$seed_thr[seeds$rep == x & seeds$sim == y]
                    
                    # Find the threshold(s)
                    thr_real_lvl_MLN_f(sim_data_k = temp_data, 
                                       p_d1_cond_ly_k = temp_p_d1_cond_ly,
                                       seed_k = temp_seed)
                  })
         })

# Save
save(thresholds_real_lvl_MLN, file = "thresholds_real_lvl_MLN.rda")

# 5 Results presented in the main text ----

# . 5.1 Threshold proportions ----

# Clear the environment
rm(list = ls())

# Import the thresholds
load("thresholds_real_lvl_MaxMI22_smooth.rda")
load("thresholds_real_lvl_CUMP.rda")
load("thresholds_real_lvl_MLN.rda")

# Extract the MaxMI-22 (smooth) thresholds
thresholds_MaxMI22 <- bind_rows(lapply(thresholds_real_lvl_MaxMI22_smooth, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")

# Extract the CUMP thresholds
thresholds_CUMP <- bind_rows(lapply(thresholds_real_lvl_CUMP, function(x) bind_rows(x, .id = "sim")), .id = "rep")

# Extract the MLN thresholds
thresholds_MLN <- bind_rows(lapply(thresholds_real_lvl_MLN, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")

# Combine the thresholds to one data frame
thresholds <- 
  thresholds_MaxMI22 %>% mutate(method = "MaxMI-22") %>%
  bind_rows(thresholds_CUMP %>% mutate(method = "CUMP")) %>%
  bind_rows(thresholds_MLN %>% mutate(method = "MLN"))

# Compute percentages of thresholds across all repetitions and conditions
thresholds %>%
  group_by(method) %>%
  summarize(N = length(na.omit(thr)), .groups = "drop") %>%
  mutate(P = 100*N/(100*80))

# Additionally, for MLN, determine what percentages of the missing cases happened because of...
# A) the two-component mixture model not fitting the data better than a one-component model,
# B) the estimated parameters not aligning with our conceptual model underlying the RTs,
# C) the estimated RT density not exhibiting bimodality.

# Extract the information:
mln_conditions <- bind_rows(lapply(thresholds_real_lvl_MLN, function(x) bind_rows(lapply(x, function(y) y$threshold_computation_conditions), .id = "sim")), .id = "rep")

# Compute percentages:
mln_conditions %>%
  # Filter out the succesfull cases
  filter(!condA | !condB | !condC) %>%
  # Create groupings:
  mutate(A = ifelse(!condA, 1, 0), # All cases where MLN failed because of A
         B = ifelse(condA & !condB, 1, 0), # All cases where MLN failed because of B
         C = ifelse(condA & condB & !condC, 1, 0)) %>% # All cases where MLN failed because of C
  group_by(A, B, C) %>%
  summarize(N = n()) %>%
  mutate(P = 100*N/(100*80)) %>%
  arrange(!A, !B, !C)

# . 5.2 Figure 5 in the main text ----

# . . 5.2.1 Setup ----

# Clear the environment
rm(list = ls())

# Import the parameters
load("params_real_lvl.rda")

# Import the conditional marginal probabilities of a correct response for SB test-takers at the population level
load("p_u1_cond_d1_real_lvl_df.rda")

# Import the thresholds
load("thresholds_real_lvl_MinMR.rda")
load("thresholds_real_lvl_MaxMI22_smooth.rda")
load("thresholds_real_lvl_CUMP.rda")
load("thresholds_real_lvl_MLN.rda")

# Extract the thresholds that minimise MR
thresholds_real_lvl_MinMR <- bind_rows(lapply(thresholds_real_lvl_MinMR, function(x) bind_rows(lapply(x, function(y) data.frame(y)), .id = "sim")), .id = "rep")
thresholds_real_lvl_MinMR <-
  thresholds_real_lvl_MinMR %>%
  rename(thr_MinMR = thr,
         mr_s_MinMR = mr_s,
         mr_p_MinMR = mr_p)

# Extract the MaxMI-22 (smooth) thresholds
thresholds_real_lvl_MaxMI22_smooth <- bind_rows(lapply(thresholds_real_lvl_MaxMI22_smooth, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")
thresholds_real_lvl_MaxMI22_smooth <-
  thresholds_real_lvl_MaxMI22_smooth %>%
  rename(thr_method = thr, 
         mr_s_method = mr_s,
         mr_p_method = mr_p)

# Extract the CUMP thresholds
thresholds_real_lvl_CUMP <- bind_rows(lapply(thresholds_real_lvl_CUMP, function(x) bind_rows(x, .id = "sim")), .id = "rep")
thresholds_real_lvl_CUMP <-
  thresholds_real_lvl_CUMP %>%
  rename(thr_method = thr, 
         mr_s_method = mr_s,
         mr_p_method = mr_p)

# Extract the MLN thresholds
thresholds_real_lvl_MLN <- bind_rows(lapply(thresholds_real_lvl_MLN, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")
thresholds_real_lvl_MLN <-
  thresholds_real_lvl_MLN %>%
  rename(thr_method = thr, 
         mr_s_method = mr_s,
         mr_p_method = mr_p)

# Combine the thresholds to one data frame
thresholds <- 
  thresholds_real_lvl_MaxMI22_smooth %>% mutate(method = "MaxMI-22") %>%
  bind_rows(thresholds_real_lvl_CUMP %>% mutate(method = "CUMP")) %>%
  bind_rows(thresholds_real_lvl_MLN %>% mutate(method = "MLN")) %>%
  full_join(thresholds_real_lvl_MinMR, by = c("rep", "sim"))

# . . 5.2.2 Draw the figure ----

# Form data frame for drawing the figure
plot_data <- 
  thresholds %>%
  mutate(rep = as.numeric(rep),
         sim = as.numeric(sim)) %>%
  # Join with parameters
  full_join(params_real_lvl %>% dplyr::select(sim, N, cor_theta_tau, iota, beta, b), 
            by = "sim") %>%
  # Join with the conditional marginal probabilities of a correct response for SB test-takers at the population level
  full_join(p_u1_cond_d1_real_lvl_df %>% mutate(sim = as.numeric(sim)), by = "sim") %>%
  # Create and compute necessary variables for drawing the figure
  mutate(# Difference between thresholds
    thr_diff = exp(thr_method) - exp(thr_MinMR),
    # Difference between MRs
    mr_diff = mr_s_method - mr_s_MinMR,
    # Signed difference between MRs
    mr_diff_sign = ifelse(thr_diff < 0, -mr_diff, mr_diff),
    # Grouping variable: the method used to compute a threshold
    method = factor(method, levels = c("MaxMI-22", "CUMP", "MLN")),
    # Grouping variable: sample size
    N_f = factor(N, labels = c("N = 500", "N = 2,500")),
    # Grouping variable: direction of association between ability and speed
    dir_assoc = ifelse(cor_theta_tau < 0,
                       "Negative association between ability and speed", 
                       "Positive association between ability and speed"),
    dir_assoc = factor(dir_assoc, levels = c("Negative association between ability and speed", 
                                             "Positive association between ability and speed")),
    # Grouping variable: RRB percentage
    rrb_perc = ifelse(iota < -4, "RRB % = 5%", "RRB % = 10%"),
    rrb_perc = factor(rrb_perc, levels = c("RRB % = 5%", "RRB % = 10%")),
    # Grouping variable: RT mode in the SB component
    sb_mode = ifelse(beta < 3.5, "SB mode = 22 sec.", "SB mode = 60 sec."),
    sb_mode = factor(sb_mode, levels = c("SB mode = 22 sec.", "SB mode = 60 sec."))) %>%
  group_by(method, dir_assoc, rrb_perc, sb_mode, N_f, p_u1_cond_d1) %>%
  summarize(Md = median(mr_diff_sign, na.rm = T), .groups = "drop")

# The plot
(p <- 
    ggplot(plot_data, mapping = aes(x = p_u1_cond_d1, y = Md, color = method, shape = method, group = method))+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = .25, linetype = "dashed")+
    geom_line(position = position_dodge(width = .05),
              linewidth = .5)+
    geom_point(position = position_dodge(width = .05),
               size = 1.5, stroke = .8, fill = "white")+
    scale_color_manual(values = c("orange", "blue", "red"))+
    scale_shape_manual(values = c(21,22,23))+
    facet_nested(rrb_perc+sb_mode~dir_assoc+N_f)+
    labs(x = "Conditional marginal probability of a correct response for SB test-takers", 
         y = expression(bar(italic(D))[MR]))+
    theme_bw()+
    theme(title = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 14)))

# Save the figure
ggsave(filename = "fig5.jpg",
       plot = p,
       width = 10000,
       height = 10000,
       units = "px",
       dpi = 1200)

# 6 Results presented in the Supplement ----

# . 6.1 Figure 2 in the Supplement ----

# . . 6.1.1 Setup ----

# Clear the environment
rm(list = ls())

# Import the parameters
load("params_real_lvl.rda")

# Import the conditional marginal probabilities of a correct response for SB test-takers at the population level
load("p_u1_cond_d1_real_lvl_df.rda")

# Import the population-level thresholds
load("thresholds_pop_lvl_MaxMI22_for_simulations.rda")

# Extract the pop-lvl thresholds
thresholds_pop_lvl_MaxMI22_for_simulations <- 
  bind_rows(lapply(thresholds_pop_lvl_MaxMI22_for_simulations, function(x) data.frame(thr_pop = x)), .id = "sim")

# Import the realised level thresholds
load("thresholds_real_lvl_MaxMI22.rda")
load("thresholds_real_lvl_MaxMI22_smooth.rda")

# Extract the MaxMI-22 thresholds
thresholds_real_lvl_MaxMI22 <- bind_rows(lapply(thresholds_real_lvl_MaxMI22, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")

# Extract the MaxMI-22 (smooth) thresholds
thresholds_real_lvl_MaxMI22_smooth <- bind_rows(lapply(thresholds_real_lvl_MaxMI22_smooth, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")

# Combine the thresholds to one data frame
thresholds <- 
  thresholds_real_lvl_MaxMI22 %>% mutate(method = "MaxMI-22 without smoothing") %>%
  bind_rows(thresholds_real_lvl_MaxMI22_smooth %>% mutate(method = "MaxMI-22 with smoothing")) %>%
  full_join(thresholds_pop_lvl_MaxMI22_for_simulations, by = "sim")

#  . . 6.1.2 Draw the figure ----

# Form data frame for drawing the figure
plot_data <- 
  thresholds %>%
  mutate(rep = as.numeric(rep),
         sim = as.numeric(sim)) %>%
  # Join with parameters
  full_join(params_real_lvl %>% dplyr::select(sim, N, cor_theta_tau, iota, beta, b), 
            by = "sim") %>%
  # Join with the conditional marginal probabilities of a correct response for SB test-takers at the population level
  full_join(p_u1_cond_d1_real_lvl_df %>% mutate(sim = as.numeric(sim)), by = "sim") %>%
  # Create and compute necessary variables for drawing the figure
  mutate(# Difference between thresholds
          thr_diff = exp(thr) - exp(thr_pop),
          # Grouping variable: the method used to compute a threshold
          method = factor(method, levels = c("MaxMI-22 without smoothing", "MaxMI-22 with smoothing")),
          # Grouping variable: sample size
          N_f = factor(N, labels = c("N = 500", "N = 2,500")),
          # Grouping variable: direction of association between ability and speed
          dir_assoc = ifelse(cor_theta_tau < 0,
                             "Negative association between ability and speed", 
                             "Positive association between ability and speed"),
          dir_assoc = factor(dir_assoc, levels = c("Negative association between ability and speed", 
                                                   "Positive association between ability and speed")),
          # Grouping variable: RRB percentage
          rrb_perc = ifelse(iota < -4, "RRB % = 5%", "RRB % = 10%"),
          rrb_perc = factor(rrb_perc, levels = c("RRB % = 5%", "RRB % = 10%")),
          # Grouping variable: RT mode in the SB component
          sb_mode = ifelse(beta < 3.5, "SB mode = 22 sec.", "SB mode = 60 sec."),
          sb_mode = factor(sb_mode, levels = c("SB mode = 22 sec.", "SB mode = 60 sec."))) %>%
  group_by(method, dir_assoc, rrb_perc, sb_mode, N_f, p_u1_cond_d1) %>%
  summarize(Md = median(thr_diff, na.rm = T), 
            Q_25 = quantile(thr_diff, probs = .25, na.rm = T),
            Q_75 = quantile(thr_diff, probs = .75, na.rm = T),
            .groups = "drop")

# The figure
(p <- 
  ggplot(plot_data,
       mapping = aes(x = p_u1_cond_d1, color = method, shape = method, group = method))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = .25, linetype = "dashed")+
  geom_line(mapping = aes(y = Md), 
            position = position_dodge(width = .02),
            linewidth = .5)+
  geom_segment(mapping = aes(y = Q_25, yend = Q_75),
               position = position_dodge(width = .02))+
  geom_point(mapping = aes(y = Md), 
             position = position_dodge(width = .02),
             size = 1.5, stroke = .8, fill = "white")+
  scale_color_manual(values = c("black", "orange"))+
  scale_shape_manual(values = c(25,21))+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  ggh4x::facet_nested(rrb_perc+sb_mode~dir_assoc+N_f)+
  labs(x = "Conditional marginal probability of a correct response for SB test-takers", 
       y = "Difference between realized and population-level threshold values")+
  theme_bw()+
  theme(title = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)))

# save the figure:
ggsave(filename = "fig2_supp.jpg",
       plot = p,
       width = 10000,
       height = 10000,
       units = "px",
       dpi = 1200)

# . 6.2 Figure 3 in the Supplement ----

# . . 6.2.1 Setup ----

# Clear the environment
rm(list = ls())

# Import the parameters
load("params_real_lvl.rda")

# Import the conditional marginal probabilities of a correct response for SB test-takers at the population level
load("p_u1_cond_d1_real_lvl_df.rda")

# Import the thresholds
load("thresholds_real_lvl_CUMP.rda")
load("thresholds_real_lvl_MLN.rda")

# Extract the CUMP thresholds
thresholds_real_lvl_CUMP <- bind_rows(lapply(thresholds_real_lvl_CUMP, function(x) bind_rows(x, .id = "sim")), .id = "rep")

# Extract the MLN thresholds
thresholds_real_lvl_MLN <- bind_rows(lapply(thresholds_real_lvl_MLN, function(x) bind_rows(lapply(x, function(y) y$thr_df), .id = "sim")), .id = "rep")

# Combine the thresholds to one data frame
thresholds <- 
  thresholds_real_lvl_CUMP %>% mutate(method = "CUMP") %>%
  bind_rows(thresholds_real_lvl_MLN %>% mutate(method = "MLN"))

# . . 6.2.2 Draw the figure ----

# Form data frame for drawing the figure
plot_data <- 
  thresholds %>%
  mutate(rep = as.numeric(rep),
         sim = as.numeric(sim)) %>%
  # Join with parameters
  full_join(params_real_lvl %>% dplyr::select(sim, N, cor_theta_tau, iota, beta, b), 
            by = "sim") %>%
  # Join with the conditional marginal probabilities of a correct response for SB test-takers at the population level
  full_join(p_u1_cond_d1_real_lvl_df %>% mutate(sim = as.numeric(sim)), by = "sim") %>%
  # Create and compute necessary variables for drawing the figure
  mutate(
    # Grouping variable: the method used to compute a threshold
    method = factor(method, levels = c("MaxMI-22", "CUMP", "MLN")),
    # Grouping variable: sample size
    N_f = factor(N, labels = c("N = 500", "N = 2,500")),
    # Grouping variable: direction of association between ability and speed
    dir_assoc = ifelse(cor_theta_tau < 0,
                       "Negative association between ability and speed", 
                       "Positive association between ability and speed"),
    dir_assoc = factor(dir_assoc, levels = c("Negative association between ability and speed", 
                                             "Positive association between ability and speed")),
    # Grouping variable: RRB percentage
    rrb_perc = ifelse(iota < -4, "RRB % = 5%", "RRB % = 10%"),
    rrb_perc = factor(rrb_perc, levels = c("RRB % = 5%", "RRB % = 10%")),
    # Grouping variable: RT mode in the SB component
    sb_mode = ifelse(beta < 3.5, "SB mode = 22 sec.", "SB mode = 60 sec."),
    sb_mode = factor(sb_mode, levels = c("SB mode = 22 sec.", "SB mode = 60 sec."))) %>%
  group_by(method, dir_assoc, rrb_perc, sb_mode, N_f, p_u1_cond_d1) %>%
  summarize(N = length(na.omit(thr)), .groups = "drop") %>%
  mutate(P = N/100)

# The plot
(p <- 
    ggplot(plot_data, mapping = aes(x = p_u1_cond_d1, y = P, color = method, shape = method, group = method))+
    geom_vline(xintercept = .25, linetype = "dashed")+
    geom_line(position = position_dodge(width = .05),
              linewidth = .5)+
    geom_point(position = position_dodge(width = .05),
               size = 1.5, stroke = .8, fill = "white")+
    scale_color_manual(values = c("blue", "red"))+
    scale_shape_manual(values = c(22,23))+
    scale_y_continuous(breaks = c(0,.25,.5,.75,1))+
    facet_nested(rrb_perc+sb_mode~dir_assoc+N_f)+
    labs(x = "Conditional marginal probability of a correct response for SB test-takers", 
         y = "Proportion of computed thresholds of all repetitions")+
    theme_bw()+
    theme(title = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 14)))

# Save the figure
ggsave(filename = "fig3_supp.jpg",
       plot = p,
       width = 10000,
       height = 10000,
       units = "px",
       dpi = 1200)

# END OF SCRIPT ----