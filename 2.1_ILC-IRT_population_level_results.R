# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for reproducing the results of employing the MaxMI-22
# method in the ILC-IRT framework at the population level

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
library(cubature)

# Settings for tidylog
crayon = function(x) cat(green(x), sep = "\n")
options("tidylog.display" = list(crayon))
rm(crayon)

# 2 Form a data frame of the person and item parameters ----

# See Table 1 in the Supplement

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
b <- sort(seq(-4, 4, length.out = 100))

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
params_pop_lvl <- expand.grid(mu_phi, mu_theta, mu_tau, 
                              var_phi, var_theta, var_tau, 
                              cor_phi_theta, cor_phi_tau, cor_theta_tau,
                              iota,
                              mu_c, sigma_c, alpha, beta,
                              g, a, b, c)

# Rename columns
colnames(params_pop_lvl) <- c("mu_phi", "mu_theta", "mu_tau",
                              "var_phi", "var_theta", "var_tau",
                              "cor_phi_theta", "cor_phi_tau", "cor_theta_tau",
                              "iota",
                              "mu_c", "sigma_c", "alpha", "beta",
                              "g", "a", "b", "c")

# Compute also the covariances between the person parameters
params_pop_lvl <-
  params_pop_lvl %>%
  mutate(cov_phi_theta = cor_phi_theta*sqrt(var_phi*var_theta),
         cov_phi_tau = cor_phi_tau*sqrt(var_phi*var_tau),
         cov_theta_tau = cor_theta_tau*sqrt(var_theta*var_tau))

# Arrange the data frame and add item number
params_pop_lvl <-
  params_pop_lvl %>%
  arrange(cor_theta_tau, iota, beta, b) %>%
  mutate(item = 1:n()) %>%
  relocate(item, .before = mu_phi)

# save
save(params_pop_lvl, file = "params_pop_lvl.rda")

# 3 Compute thresholds for RRB identification ----

# . 3.1 Compute thresholds based on minimising misclassification rate of X ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_pop_lvl.rda")

# Read the user-defined functions from the "2.0_ILC-IRT_population_level_functions.R" script
source("2.0_ILC-IRT_population_level_functions.R")

# Note! Here we do not have go through all 800 items in data frame "params_pop_lvl".
# We see this from the formula for MR (see Eq. 30 in the Supplement): 
# MR = P(Delta = 0, X = 0 | iota, mu_c, sigma_c) + P(Delta = 1, X = 1 | iota, alpha, beta)
# That is, the only item parameters that affect MR are iota, mu_c, sigma_c, alpha, and beta
# However, out of these, only iota and beta were varied.

# First, extract distinct conditions from "params_pop_lvl" based on iota and beta:
temp <- 
  params_pop_lvl %>% 
  distinct(iota, beta, .keep_all = T) %>%
  select(item, mu_phi:var_tau, cov_phi_theta:cov_theta_tau, iota, beta, alpha, mu_c, sigma_c)
# temp should have 4 rows (2 unique values of iota and 2 unique values of beta, i.e., 2x2=4 distinct conditions)

# Second, apply "thr_pop_lvl_MinMR_f" over the four distinct conditions
thresholds_pop_lvl_MinMR <- lapply(temp$item, function(x) thr_pop_lvl_MinMR_f(k = x, params = temp))

# Save the thresholds
save(thresholds_pop_lvl_MinMR, file = "thresholds_pop_lvl_MinMR.rda")

# . 3.2 MaxMI-22 ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_pop_lvl.rda")

# Read the user-defined functions from the "2.0_ILC-IRT_population_level_functions.R" script
source("2.0_ILC-IRT_population_level_functions.R")

# Apply "thr_pop_lvl_MaxMI22_f" over the 800 different conditions
# This may take a while so you may want to do it in item blocks
thresholds_pop_lvl_MaxMI22 <- lapply(params_pop_lvl$item, function(x) thr_pop_lvl_MaxMI22_f(k = x, params = params_pop_lvl))

# When finished, combine to a data frame
thresholds_pop_lvl_MaxMI22 <- bind_rows(thresholds_pop_lvl_MaxMI22, .id = "item")

# Save the thresholds
save(thresholds_pop_lvl_MaxMI22, file = "thresholds_pop_lvl_MaxMI22.rda")

# 4 Compute the conditional marginal probabilities of a correct response for SB test-takers ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_pop_lvl.rda")

# Read the user-defined functions from the "2.0_ILC-IRT_population_level_functions.R" script
source("2.0_ILC-IRT_population_level_functions.R")

# Apply "p_u1_cond_d1_f" over the 800 different conditions
res <- lapply(params_pop_lvl$item, function(x) p_u1_cond_d1_f(k = x, params = params_pop_lvl))

# When finished, combine to a data frame
p_u1_cond_d1_pop_lvl_df <- bind_rows(lapply(res, function(x) data.frame(p_u1_cond_d1 = x)), .id = "item")

# save
save(p_u1_cond_d1_pop_lvl_df, file = "p_u1_cond_d1_pop_lvl_df.rda")


# 5 Draw Figure 4 in the main text ----

# . 5.1 Setup ----

# Clear the environment
rm(list = ls())

# Import the parameter data frame
load("params_pop_lvl.rda")

# Import the thresholds
load("thresholds_pop_lvl_MinMR.rda")
load("thresholds_pop_lvl_MaxMI22.rda")

# Import the conditional marginal probabilities
load("p_u1_cond_d1_pop_lvl_df.rda")

# Combine the thresholds obtained by minimising the MR to a data frame
# Also, expand the data frame so that each of the 800 conditions gets a threshold
# Recall that the only parameters affecting these thresholds were iota and beta.
# Conditions 1-100: iota = -4.27 and beta = 3.25
# Conditions 101-200: iota = -3.29 and beta = 3.25
# Conditions 201-300: iota = -4.27 and beta = 4.25
# Conditions 301-400: iota = -3.29 and beta = 4.25
# Then the same repeats for conditions 401-800
thresholds_pop_lvl_MinMR <-
  data.frame(item = 1:800,
             thr_MinMR = c(rep(thresholds_pop_lvl_MinMR[[1]]$thr, 100), 
                           rep(thresholds_pop_lvl_MinMR[[2]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[3]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[4]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[1]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[2]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[3]]$thr, 100),
                           rep(thresholds_pop_lvl_MinMR[[4]]$thr, 100)),
             mr_MinMR = c(rep(thresholds_pop_lvl_MinMR[[1]]$mr, 100), 
                          rep(thresholds_pop_lvl_MinMR[[2]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[3]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[4]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[1]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[2]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[3]]$mr, 100),
                          rep(thresholds_pop_lvl_MinMR[[4]]$mr, 100)))

# Mutate the column names in "thresholds_pop_lvl_MaxMI22"
thresholds_pop_lvl_MaxMI22 <- thresholds_pop_lvl_MaxMI22 %>% rename(thr_MaxMI22 = thr, mr_MaxMI22 = mr)

# Combine the thresholds to one data frame
thresholds <- 
  thresholds_pop_lvl_MaxMI22 %>%
  full_join(thresholds_pop_lvl_MinMR, by = "item")

# Transform the thresholds to the "original" scale t = exp(l)
thresholds <-
  thresholds %>%
  mutate(thr_MaxMI22 = exp(thr_MaxMI22),
         thr_MinMR = exp(thr_MinMR))
  
# . 5.2 Draw the figure ----

# Form the data frame for drawing the figure
plot_data <-
  thresholds %>%
  # Add beta, iota, and cor_theta_tau
  full_join(params_pop_lvl %>% select(item, beta, iota, cor_theta_tau), by = "item") %>%
  # Add the conditional marginal probabilities
  full_join(p_u1_cond_d1_pop_lvl_df, by = "item") %>%
  # Create and compute necessary variables for drawing the figure
  mutate(# Difference between thresholds
         diff_t = thr_MaxMI22 - thr_MinMR,
         # Difference between MRs
         diff_mr = mr_MaxMI22 - mr_MinMR,
         # Signed difference between MRs:
         diff_mr_sign = ifelse(diff_t < 0, -diff_mr, diff_mr),
         # Grouping variable: direction of association between ability and speed
         dir_assoc = ifelse(cor_theta_tau < 0,
                            "Negative association between ability and speed", 
                            "Positive association between ability and speed"),
         dir_assoc = factor(dir_assoc, levels = c("Negative association between ability and speed", 
                                                  "Positive association between ability and speed")),
         # Grouping variable: distinct combinations of iota and beta
         iota_beta = ifelse(iota < -4 & beta < 3.5,
                            "RRB % = 5%; SB mode = 22 sec.",
                            NA),
         iota_beta = ifelse(iota < -4 & beta > 4,
                            "RRB % = 5%; SB mode = 60 sec.",
                            iota_beta),
         iota_beta = ifelse(iota > -4 & beta < 3.5,
                            "RRB % = 10%; SB mode = 22 sec.",
                            iota_beta),
         iota_beta = ifelse(iota > -4 & beta > 4,
                            "RRB % = 10%; SB mode = 60 sec.",
                            iota_beta),
         iota_beta = factor(iota_beta,
                            levels = c("RRB % = 5%; SB mode = 22 sec.",
                                       "RRB % = 5%; SB mode = 60 sec.",
                                       "RRB % = 10%; SB mode = 22 sec.",
                                       "RRB % = 10%; SB mode = 60 sec.")))


# Draw the figure:

# Plot A
(p_A <-
    ggplot()+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = .25, linetype = "dashed")+
    geom_line(plot_data %>% filter(p_u1_cond_d1 < .25), 
              mapping = aes(x = p_u1_cond_d1, y = diff_t, 
                            color = iota_beta, linetype = iota_beta),
              linewidth = .8)+
    geom_line(plot_data %>% filter(p_u1_cond_d1 > .25), 
              mapping = aes(x = p_u1_cond_d1, y = diff_t, 
                            color = iota_beta, linetype = iota_beta),
              linewidth = .8)+
    facet_wrap(~dir_assoc)+
    scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00"))+
    scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00"))+
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))+
    scale_shape_manual(values = c(21, 22, 23, 24))+
    labs(x = "Conditional marginal probability of a correct response for SB test-takers", 
         y = expression(paste(italic(D)[italic(T)], " (sec.)")),
         title = "A")+
    theme_bw()+
    theme(title = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.key.width = unit(3, "line"),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE),
           color=guide_legend(nrow=2,byrow=TRUE)))

# Plot B
(p_B <-
    ggplot()+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = .25, linetype = "dashed")+
    geom_line(plot_data %>% filter(p_u1_cond_d1 < .25), 
              mapping = aes(x = p_u1_cond_d1, y = diff_mr_sign, 
                            color = iota_beta, linetype = iota_beta),
              linewidth = .8)+
    geom_line(plot_data %>% filter(p_u1_cond_d1 > .25), 
              mapping = aes(x = p_u1_cond_d1, y = diff_mr_sign, 
                            color = iota_beta, linetype = iota_beta),
              linewidth = .8)+
    facet_wrap(~dir_assoc)+
    scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00"))+
    scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00"))+
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))+
    scale_shape_manual(values = c(21, 22, 23, 24))+
    labs(x = "Conditional marginal probability of a correct response for SB test-takers", 
         y = expression(italic(D)[MR]),
         title = "B")+
    theme_bw()+
    theme(title = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.key.width = unit(3, "line"),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE),
           color=guide_legend(nrow=2,byrow=TRUE)))

# The final plot
(p <- ggarrange(p_A, p_B, nrow = 2, legend = "bottom", common.legend = T))

# Save the figure
ggsave(filename = "fig4.jpg",
       plot = p,
       width = 10000,
       height = 10000,
       units = "px",
       dpi = 1200)

# END OF SCRIPT ----