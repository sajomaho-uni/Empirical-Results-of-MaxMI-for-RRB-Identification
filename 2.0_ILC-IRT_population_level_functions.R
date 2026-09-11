# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for user-defined functions that are used in the 
# "2.1_ILC-IRT_population_level_results.R" script.

# Author information:
# Santeri Holopainen
# Turku Research Institute for Learning Analytics (TRILA)
# University of Turku

# Date: 11.9.2026

# R version: 4.5.0

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# Function argument explanations:

# Person parameters:
# phi = latent engagement tendency of a person
# theta = latent ability of a person
# tau = latent speed of a person
# person_par = a 2x1 vector (phi, theta)
# mu_P = a 3x1 mean vector of the person parameters (phi, theta, and tau)
# Sigma_P = a 3x3 covariance matrix of the person parameters (phi, theta, and tau)

# Item parameters:
# iota = item engagement difficulty parameter
# a = item discrimination parameter
# b = item difficulty parameter
# c = item pseudo-guessing parameter
# alpha = time discrimination parameter of the item
# beta = time intensity parameter of the item
# g = chance level of the item
# mu_c = common mean for disengaged log RTs
# sigma_c = common standard deviation for disengaged log RTs

# Other arguments:
# l0 = the log-time threshold
# eps = a very small numeric value (eps = .Machine$double.xmin by default)
# p_d1 = P(Delta = 1 | iota)
# p_u1 = P(U = 1 | iota, a, b, c, g) = P(U = 1, Delta = 1 | iota, a, b, c) + P(U = 1, Delta = 0 | g)
# k = item identifier
# params = data frame of person and item parameters

# 1 IRT models for engagement and response probabilities ----

# Probability that person answers under engaged state (Rasch, i.e., 1PL IRT)
p_engagement <- function(phi, iota){
  plogis(phi-iota)
}

# Probability of a correct response in the solution behaviour component (3PL IRT)
p_correct <- function(theta, a, b, c){
  # If a = 1 and c = 0, this is the Rasch model
  c + (1-c)*plogis(a*(theta-b))
}

# 2 Conditional means and standard deviations of the person parameters ----

# These formulas assume that the person parameters are multivariate normally distributed.
# See Eqs. 11 and 12 in the Supplement.
# mu_P and Sigma_P shall be formed before using the following functions.

# E(tau | phi)
mu_tau_cond_phi <- function(phi, mu_P, Sigma_P){
  mu_P[3] + (Sigma_P[3,1] / Sigma_P[1,1]) * (phi - mu_P[1])
}

# SD(tau | phi)
sigma_tau_cond_phi <- function(Sigma_P){
  sqrt(Sigma_P[3,3] - Sigma_P[3,1]^2 / Sigma_P[1,1])
}

# E(theta | phi)
mu_theta_cond_phi <- function(phi, mu_P, Sigma_P){
  mu_P[2] + (Sigma_P[2,1] / Sigma_P[1,1]) * (phi - mu_P[1])
}

# SD(theta | phi)
sigma_theta_cond_phi <- function(Sigma_P){
  sqrt(Sigma_P[2,2] - Sigma_P[2,1]^2 / Sigma_P[1,1])
}

# E(tau | phi, theta)
mu_tau_cond_phi_theta <- function(phi, theta, mu_P, Sigma_P){
  mu_P[3] + Sigma_P[3,1:2] %*% solve(Sigma_P[1:2,1:2]) %*% (c(phi, theta) - mu_P[1:2])
}

# SD(tau | phi, theta)
sigma_tau_cond_phi_theta <- function(Sigma_P){
  sqrt(Sigma_P[3,3] - Sigma_P[3,1:2] %*% solve(Sigma_P[1:2,1:2]) %*% Sigma_P[1:2,3])
}

# 3 Functions to be integrated numerically ----

# In this section, we define functions that shall be integrated numerically when
# acquiring the population-level results of employing the MaxMI-22 method.

# In the following integrands, the person parameters are the variables of integration.
# See the Supplement for the formulas.

# 1) Integrand for computing the marginal probability P(Delta = 1 | iota).
# See Eq. 13 in the Supplement.
int_p_d1 <- function(phi, iota, mu_P, Sigma_P){
  p_engagement(phi = phi, iota = iota)*dnorm(phi, mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))
}

# 2) Integrand for computing the marginal probability P(U = 1, Delta = 1 | iota, a, b, c)
# See Eq. 14 in the Supplement
int_p_u1d1 <- function(person_par, iota, a, b, c, mu_P, Sigma_P){
  p_engagement(phi = person_par[1], iota = iota)*
    p_correct(theta = person_par[2], a = a, b = b, c = c)*
    dnorm(person_par[1], mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
    dnorm(person_par[2], 
          mean = mu_theta_cond_phi(phi = person_par[1], mu_P = mu_P, Sigma_P = Sigma_P),
          sd = sigma_theta_cond_phi(Sigma_P = Sigma_P))
}

# 3) Integrand for computing the log RT density at l0 in the SB component
# See Eqs. 17-19 in the Supplement.
int_f_L_d1 <- function(phi, l0, iota, alpha, beta, mu_P, Sigma_P){
  p_engagement(phi = phi, iota = iota)*
    dnorm(phi, mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
    dnorm(beta - l0, 
          mean = mu_tau_cond_phi(phi = phi, mu_P = mu_P, Sigma_P = Sigma_P),
          sd = sqrt(alpha^2 + sigma_tau_cond_phi(Sigma_P = Sigma_P)^2))
}

# 4) Integrand for computing the marginal probability P(X = 1, Delta = 1 | iota, alpha, beta) at l0
# See Eqs. 20-22 in the Supplement.
int_p_x1d1 <- function(phi, l0, iota, alpha, beta, mu_P, Sigma_P){
  p_engagement(phi = phi, iota = iota)*
    dnorm(phi, mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
    pnorm((mu_tau_cond_phi(phi = phi, mu_P = mu_P, Sigma_P = Sigma_P) + l0 - beta)/sqrt(alpha^2+sigma_tau_cond_phi(Sigma_P = Sigma_P)^2))
}

# 5) Integrand for computing the marginal probability P(X = 1, U = 1, Delta = 1 | iota, a, b, c, alpha, beta) at l0
# See Eq. 27 in the Supplement.
int_p_x1u1d1 <- 
  function(person_par, l0, iota, a, b, c, alpha, beta, mu_P, Sigma_P){
    p_engagement(phi = person_par[1], iota = iota)*
      dnorm(person_par[1], mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
      p_correct(theta = person_par[2], a = a, b = b, c = c)*
      dnorm(person_par[2], 
            mean = mu_theta_cond_phi(phi = person_par[1], mu_P = mu_P, Sigma_P = Sigma_P), 
            sd = sigma_theta_cond_phi(Sigma_P = Sigma_P))*
      pnorm((mu_tau_cond_phi_theta(phi = person_par[1], theta = person_par[2], mu_P = mu_P, Sigma_P = Sigma_P) + l0 - beta)/sqrt(alpha^2 + sigma_tau_cond_phi_theta(Sigma_P = Sigma_P)^2))
  }

# 4 Functions for computing MIRT and the misclassification rate of X at l0 ----

# In the following functions, arguments p_d1 and p_u1 must be computed before using them.
# p_d1 is computed by integrating the integrand int_p_d1 to receive P(Delta = 1 | iota).
# p_u1 is computed by...
# 1) integrating the integrand int_p_u1d1 to receive P(U = 1, Delta = 1 | iota, a, b, c),
# 2) computing (1-p_d1)*g to receive P(U = 1, Delta = 0 | iota, g), and
# 3) adding P(U = 1, Delta = 1 | iota, a, b, c) and P(U = 1, Delta = 0 | iota, g) together.

# Function for computing MIRT at a specific log-time threshold l0
mirt_f <- function(l0, 
                   p_d1, p_u1, 
                   iota, 
                   mu_c, sigma_c, alpha, beta, 
                   a, b, c, g, mu_P, Sigma_P, 
                   eps = .Machine$double.xmin){
  # P(X=1)
  p_x1d1 <- integrate(int_p_x1d1, 
                      lower = -100, upper = 100,
                      l0 = l0,
                      iota = iota,
                      alpha = alpha, beta = beta,
                      mu_P = mu_P, Sigma_P = Sigma_P)$value
  p_x1d0 <- (1-p_d1)*pnorm((l0-mu_c)/sigma_c)
  p_x1 <- p_x1d1 + p_x1d0
  
  # P(X=0)
  p_x0 <- 1-p_x1
  
  # P(U=0)
  p_u0 <- 1-p_u1
  
  # P(X=1,U=1)
  p_x1u1d1 <- hcubature(int_p_x1u1d1,
                        lowerLimit = c(-100, -100),
                        upperLimit = c(100, 100),
                        l0 = l0,
                        iota = iota,
                        a = a, b = b, c = c,
                        alpha = alpha, beta = beta,
                        mu_P = mu_P, Sigma_P = Sigma_P)$integral
  p_x1u1d0 <- (1-p_d1)*pnorm((l0-mu_c)/sigma_c)*g
  p_x1u1 <- p_x1u1d1 + p_x1u1d0
  
  # P(X=1,U=0)
  p_x1u0 <- p_x1 - p_x1u1
  
  # P(X=0,U=1)
  p_x0u1 <- p_u1 - p_x1u1
  
  # P(X=0,U=0)
  p_x0u0 <- 1-p_x1-p_u1+p_x1u1
  
  # MIRT
  p_x1u1 * (log2(p_x1u1 + eps) - log2(p_x1 * p_u1 + eps)) +
    p_x1u0 * (log2(p_x1u0 + eps) - log2(p_x1 * p_u0 + eps)) +
    p_x0u1 * (log2(p_x0u1 + eps) - log2(p_x0 * p_u1 + eps)) +
    p_x0u0 * (log2(p_x0u0 + eps) - log2(p_x0 * p_u0 + eps))
}

# Function for computing the misclassification rate of X at a specific log-time threshold l0.
# MR = P(Delta = 0, X = 0 | iota, mu_c, sigma_c) + P(Delta = 1, X = 1 | iota, alpha, beta)
# See Eq. 31 in the Supplement.
MR_f <- function(l0, p_d1, mu_c, sigma_c, iota, alpha, beta, mu_P, Sigma_P){
  (1-p_d1)*(1-pnorm((l0-mu_c)/sigma_c))+
    integrate(int_p_x1d1, 
              lower = -100, upper = 100,
              l0 = l0,
              iota = iota,
              alpha = alpha, beta = beta,
              mu_P = mu_P, Sigma_P = Sigma_P)$value
}

# 5 Functions for population-level computations ----

# Function for computing a threshold by minimising the misclassification rate of X
thr_pop_lvl_MinMR_f <- function(k, params){
  
  # Filter out other items from the parameter data frame
  fun_par <- params[params$item == k,]
  
  # Extract parameters from "fun_par"
  mu_P <- c(fun_par$mu_phi,fun_par$mu_theta,fun_par$mu_tau)
  Sigma_P <- matrix(c(fun_par$var_phi,       fun_par$cov_phi_theta, fun_par$cov_phi_tau,
                      fun_par$cov_phi_theta, fun_par$var_theta,     fun_par$cov_theta_tau,
                      fun_par$cov_phi_tau,   fun_par$cov_theta_tau, fun_par$var_tau), 
                    nrow = 3, byrow = T)
  iota <- fun_par$iota
  beta <- fun_par$beta
  alpha <- fun_par$alpha
  mu_c <- fun_par$mu_c
  sigma_c <- fun_par$sigma_c
  
  # Compute P(Delta = 1) 
  p_d1 <- integrate(int_p_d1, 
                    lower = -100, upper = 100, 
                    iota = iota, 
                    mu_P = mu_P, Sigma_P = Sigma_P)$value
  
  # Compute the threshold
  thr <- 
    optimise(MR_f, 
             interval = c(0,beta - mu_P[3]),
             p_d1 = p_d1,
             mu_c = mu_c, sigma_c = sigma_c,
             iota = iota,
             beta = beta, alpha = alpha,
             mu_P = mu_P, Sigma_P = Sigma_P)$minimum
  
  # Also, compute the MR at the threshold
  mr <- MR_f(l0 = thr,
             p_d1 = p_d1,
             mu_c = mu_c, sigma_c = sigma_c,
             iota = iota,
             beta = beta, alpha = alpha,
             mu_P = mu_P, Sigma_P = Sigma_P)
  
  # Done
  list(thr = thr, mr = mr)
}

# Function for computing a MaxMI-22 threshold
thr_pop_lvl_MaxMI22_f <- function(k, params){
  
  # Filter out other items from the parameter data frame
  fun_par <- params[params$item == k,]
  
  # Extract parameters from "fun_par"
  mu_P <- c(fun_par$mu_phi,fun_par$mu_theta,fun_par$mu_tau)
  Sigma_P <- matrix(c(fun_par$var_phi,       fun_par$cov_phi_theta, fun_par$cov_phi_tau,
                      fun_par$cov_phi_theta, fun_par$var_theta,     fun_par$cov_theta_tau,
                      fun_par$cov_phi_tau,   fun_par$cov_theta_tau, fun_par$var_tau), 
                    nrow = 3, byrow = T)
  iota <- fun_par$iota
  a <- fun_par$a
  b <- fun_par$b
  c <- fun_par$c
  g <- fun_par$g
  beta <- fun_par$beta
  alpha <- fun_par$alpha
  mu_c <- fun_par$mu_c
  sigma_c <- fun_par$sigma_c
  
  # Compute P(Delta = 1) 
  p_d1 <- integrate(int_p_d1, 
                    lower = -100, upper = 100, 
                    iota = iota, 
                    mu_P = mu_P, Sigma_P = Sigma_P)$value
  
  # Compute P(Y = 1)
  p_u1 <- (1-p_d1)*g + hcubature(int_p_u1d1,
                                 lowerLimit = c(-Inf, -Inf),
                                 upperLimit = c(Inf, Inf),
                                 iota = iota,
                                 a = a, b = b, c = c,
                                 mu_P = mu_P, Sigma_P = Sigma_P)$integral
  
  # Before computing the threshold, we restrict the search space of the global maximum of MIRT
  # We do this by conducting a "rough" iterative search for a log-time point (l2)
  # where MIRT is smaller than MIRT at an earlier log-time point (l1)
  
  # Starting points:
  l1 <- 1; l2 <- 1.1
  
  # MIRT at l1
  mirt_at_l1 <- mirt_f(l0 = l1,
                       p_d1 = p_d1, p_u1 = p_u1,
                       iota = iota,
                       mu_c = mu_c, sigma_c = sigma_c,
                       beta = beta, alpha = alpha,
                       a = a, b = b, c = c, g = g,
                       mu_P = mu_P, Sigma_P = Sigma_P)
  
  # MIRT at l2
  mirt_at_l2 <- mirt_f(l0 = l2,
                       p_d1 = p_d1, p_u1 = p_u1,
                       iota = iota,
                       mu_c = mu_c, sigma_c = sigma_c,
                       beta = beta, alpha = alpha,
                       a = a, b = b, c = c, g = g,
                       mu_P = mu_P, Sigma_P = Sigma_P)
  
  # Iterations with a while loop.
  # Repeat until MIRT at l2 is less than or equal to MIRT at l1:
  while(mirt_at_l2 > mirt_at_l1){
    # Update the time points.
    # l1 becomes l2 and l2 is increased by adding 0.1 to it:
    l1 <- l2; l2 <- l2 + 0.1
    
    # Update the MIRT values.
    # MIRT at l1 becomes MIRT at l2.
    # MIRT at l2 is computed again:
    mirt_at_l1 <- mirt_at_l2
    mirt_at_l2 <- mirt_f(l0 = l2,
                         p_d1 = p_d1, p_u1 = p_u1,
                         iota = iota,
                         mu_c = mu_c, sigma_c = sigma_c,
                         beta = beta, alpha = alpha,
                         a = a, b = b, c = c, g = g,
                         mu_P = mu_P, Sigma_P = Sigma_P)
    
  }
  # After the loop is done, we should have a suitable upper boundary (l2)
  # for the search space of the global maximum of MIRT
  
  # Compute the MaxMI-22 threshold
  thr <- 
    optimise(mirt_f, 
             interval = c(1, l2),
             p_d1 = p_d1, p_u1 = p_u1,
             iota = iota,
             mu_c = mu_c, sigma_c = sigma_c,
             beta = beta, alpha = alpha,
             a = a, b = b, c = c, g = g,
             mu_P = mu_P, Sigma_P = Sigma_P,
             maximum = T)$maximum
  
  # Also, compute the MR at the threshold
  mr <- MR_f(l0 = thr,
             p_d1 = p_d1,
             mu_c = mu_c, sigma_c = sigma_c,
             iota = iota,
             beta = beta, alpha = alpha,
             mu_P = mu_P, Sigma_P = Sigma_P)
  
  # Done
  list(thr = thr, mr = mr)
}

# Function for computing the conditional marginal probability of a correct response
# for SB test-takers
p_u1_cond_d1_f <- function(k, params){
  
  # Filter out other items from the parameter data frame
  fun_par <- params[params$item == k,]
  
  # Extract parameters from "fun_par"
  mu_P <- c(fun_par$mu_phi,fun_par$mu_theta,fun_par$mu_tau)
  Sigma_P <- matrix(c(fun_par$var_phi,       fun_par$cov_phi_theta, fun_par$cov_phi_tau,
                      fun_par$cov_phi_theta, fun_par$var_theta,     fun_par$cov_theta_tau,
                      fun_par$cov_phi_tau,   fun_par$cov_theta_tau, fun_par$var_tau), 
                    nrow = 3, byrow = T)
  iota <- fun_par$iota
  a <- fun_par$a
  b <- fun_par$b
  c <- fun_par$c
  
  # Compute P(Delta = 1) 
  p_d1 <- integrate(int_p_d1, 
                    lower = -100, upper = 100, 
                    iota = iota, 
                    mu_P = mu_P, Sigma_P = Sigma_P)$value
  
  # Compute P(U = 1, Delta = 1)
  p_u1d1 <- hcubature(int_p_u1d1,
                      lowerLimit = c(-100, -100),
                      upperLimit = c(100, 100),
                      iota = iota,
                      a = a, b = b, c = c,
                      mu_P = mu_P, Sigma_P = Sigma_P)$integral
  
  # Compute P(U = 1 | Delta = 1)
  p_u1_cond_d1 <- p_u1d1 / p_d1
  
  # Done
  p_u1_cond_d1
}

# END OF SCRIPT ----