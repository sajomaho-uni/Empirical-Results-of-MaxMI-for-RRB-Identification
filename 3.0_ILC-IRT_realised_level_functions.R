# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for user-defined functions that are used in the 
# "3.1_ILC-IRT_realised_level_results.R" script.

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
# mu_phi = mean of phi
# sigma_phi = standard deviation of phi
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
# loglik = log-likelihood of the fitted model
# q = number of parameters in the model
# n = sample size
# z = the variable of interest (e.g., empirical MIRT)
# p = probability of success
# eps = a very small numeric value (eps = .Machine$double.xmin by default)
# u = the response correctness vector
# l = the log RT vector
# l0 = the log-time threshold vector (essentially, l0 = l)
# w = a vector of weights
# k = simulation condition indicator
# params = a data frame of simulation parameters
# seed = seed number for reproducing the same result
# n_gh = number of nodes in Gauss-Hermite
# sim_data = a nested list of all generated data sets
# sim_data_k = one generated data set
# p_d1_df = a pre-computed data frame of probabilities P(Delta = 1 | iota)
# p_d1_cond_ly_k = a pre-computed vector of posterior probabilities P(Delta = 1 | l, y)
# n_cut = number of observations with the smallest response times to cut from the data prior to threshold computation
### (In some cases, it is beneficial to clean the data a little bit from spurious noise in the empirical MIRT with small RTs)
# limit = should the search space for the maximum of the empirical MIRT be limited?

# 1 General functions ----

# Function for calculating a fitted model's Bayesian Information Criterion (BIC)
My_BIC <- function(loglik, q, n){
  q*log(n)-2*loglik
}

# Function for finding local minima
local_minima <- function(z){
  s <- sign(diff(z))
  mask <- s != 0
  k <- which(diff(s[mask]) == 2)
  i <- which(s != 0)
  (i[k] + i[k+1] + 1) / 2
}

# Function for computing the entropy of a binary variable with probability p
h_f <- function(p, eps = .Machine$double.xmin){
  -p*log2(pmax(p,eps))-(1-p)*log2(pmax(1-p,eps))
}

# Function for computing the empirical CUMP curve, 
# i.e., the empirical conditional probability of a correct response given X = 1
emp_cump_f <- function(u){ 
  cumsum(u) / 1:length(u)
}

# Function for computing empirical MIRT when both X and Y=U are binary
# This function returns a vector of the empirical MIRT values
emp_mirt_x2y2_f <- function(u, eps = .Machine$double.xmin){
  n <- length(u) # sample size
  p_u1 <- mean(u) # P(U = 1)
  p_x1 <- 1:n/n # P(X = 1) = P(T <= T_k)
  p_x1u1 <- (cumsum(u)/n) # P(X = 1, U = 1)
  p_x0u1 <- p_u1-p_x1u1 # P(X = 0, U = 1)
  p_u1_cond_x1 <- p_x1u1 / p_x1  # P(U = 1 | X = 1) 
  p_u1_cond_x0 <- p_x0u1 / pmax((1-p_x1), eps) # P(U = 1 | X = 0)
  h_f(p_u1)-(p_x1*h_f(p_u1_cond_x1) + (1-p_x1)*h_f(p_u1_cond_x0)) # empirical MIRT
}

# Probability that person answers under engaged state (Rasch, i.e., 1PL IRT)
p_engagement <- function(phi, iota){
  plogis(phi-iota)
}

# Probability of a correct response in the solution behaviour component (3PL IRT)
p_correct <- function(theta, a, b, c){
  # If a = 1 and c = 0, this is the Rasch model
  c + (1-c)*plogis(a*(theta-b))
}

# Integrand for pre-computing the marginal probability P(Delta = 1 | iota)
# See Eq. 12 in the Supplement.
int_p_d1 <- function(phi, iota, mu_P, Sigma_P){
  p_engagement(phi = phi, iota = iota)*dnorm(phi, mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))
}

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

# Integrand for computing the marginal probability P(U = 1, Delta = 1 | iota, a, b, c)
# See Eq. 13 in the Supplement
int_p_u1d1 <- function(person_par, iota, a, b, c, mu_P, Sigma_P){
  p_engagement(phi = person_par[1], iota = iota)*
    p_correct(theta = person_par[2], a = a, b = b, c = c)*
    dnorm(person_par[1], mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
    dnorm(person_par[2], 
          mean = mu_theta_cond_phi(phi = person_par[1], mu_P = mu_P, Sigma_P = Sigma_P),
          sd = sigma_theta_cond_phi(Sigma_P = Sigma_P))
}

# Function for computing the conditional marginal probability of a correct response
# for SB test-takers for the kth item at the population level
p_u1_cond_d1_f <- function(k, params){
  
  # Filter out other items from the parameter data frame
  fun_par <- params[params$sim == k,]
  
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

# Integrand for computing the marginal probability P(X = 1, Delta = 1 | iota, alpha, beta) at l0
# See Eqs. 19-21 in the Supplement.
int_p_x1d1 <- function(phi, l0, iota, alpha, beta, mu_P, Sigma_P){
  p_engagement(phi = phi, iota = iota)*
    dnorm(phi, mean = mu_P[1], sd = sqrt(Sigma_P[1,1]))*
    pnorm((mu_tau_cond_phi(phi = phi, mu_P = mu_P, Sigma_P = Sigma_P) + l0 - beta)/sqrt(alpha^2+sigma_tau_cond_phi(Sigma_P = Sigma_P)^2))
}

# Integrand for computing the marginal probability P(X = 1, U = 1, Delta = 1 | iota, a, b, c, alpha, beta) at l0
# See Eq. 26 in the Supplement.
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

# Function for computing MIRT at a specific log-time threshold l0 at the population level
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

# 2 Functions for unimodal monotone regression ----

# We utilised OpenAI's GPT-5.2 model to write the following functions

# Block log-likelihood for Bernoulli in a block with totals (sumwy, sumw)
block_ll_binom <- function(sumwy, sumw, eps = .Machine$double.xmin) {
  p <- sumwy / sumw
  p <- pmin(pmax(p, eps), 1 - eps)
  sumwy * log(p) + (sumw - sumwy) * log1p(-p)
}

# Function for computing isotonic (nondecreasing) MLE logLik at m for every prefix 1:m
prefix_iso_ll_binom <- function(y, w = NULL, eps = .Machine$double.xmin) {
  n <- length(y)
  if (is.null(w)) w <- rep(1, n)
  
  # stacks for blocks
  sw  <- numeric(n)
  swy <- numeric(n)
  lvl <- numeric(n)  # not strictly needed, but keeps pooling condition simple
  llb <- numeric(n)  # block loglik
  mblk <- 0L
  
  # total loglik after fitting prefix 1:i
  prefLL <- numeric(n)
  totalLL <- 0
  
  for (i in seq_len(n)) {
    # create new block (i:i)
    mblk <- mblk + 1L
    sw[mblk]  <- w[i]
    swy[mblk] <- w[i] * y[i]
    lvl[mblk] <- swy[mblk] / sw[mblk]
    llb[mblk] <- block_ll_binom(swy[mblk], sw[mblk], eps)
    
    totalLL <- totalLL + llb[mblk]
    
    # pool adjacent violators
    while (mblk >= 2L && lvl[mblk - 1L] > lvl[mblk]) {
      # remove their old contributions
      totalLL <- totalLL - llb[mblk - 1L] - llb[mblk]
      
      # merge into (mblk-1)
      sw[mblk - 1L]  <- sw[mblk - 1L]  + sw[mblk]
      swy[mblk - 1L] <- swy[mblk - 1L] + swy[mblk]
      lvl[mblk - 1L] <- swy[mblk - 1L] / sw[mblk - 1L]
      llb[mblk - 1L] <- block_ll_binom(swy[mblk - 1L], sw[mblk - 1L], eps)
      
      # add merged contribution
      totalLL <- totalLL + llb[mblk - 1L]
      
      # pop
      mblk <- mblk - 1L
    }
    
    prefLL[i] <- totalLL
  }
  
  prefLL
}

# Function for computing the isotonic (nonincreasing) MLE logLik on that suffix for every suffix i:n
suffix_antiiso_ll_binom <- function(y, w = NULL, eps = .Machine$double.xmin) {
  n <- length(y)
  if (is.null(w)) w <- rep(1, n)
  
  # Nonincreasing on i:n == nondecreasing on reversed prefix (1:(n-i+1))
  pref_rev <- prefix_iso_ll_binom(rev(y), rev(w), eps = eps)
  
  # sufLL[i] = loglik of fit on i:n
  # corresponds to pref_rev[n - i + 1]
  sufLL <- numeric(n)
  for (i in seq_len(n)) sufLL[i] <- pref_rev[n - i + 1L]
  sufLL
}


# Function for unimodal fit: choose mode location by maximum likelihood.
unimodal_regression <- function(L, y, w = NULL, eps = .Machine$double.xmin) {
  n <- length(y)
  if (length(L) != n) stop("L and y must have same length")
  if (is.unsorted(L, strictly = FALSE)) stop("L must be sorted ascending")
  if (!all(y %in% c(0,1))) stop("y must be 0/1")
  if (is.null(w)) w <- rep(1, n)
  if (length(w) != n) stop("w must have same length as y")
  
  prefLL <- prefix_iso_ll_binom(y, w, eps = eps)
  sufLL  <- suffix_antiiso_ll_binom(y, w, eps = eps)
  totalLL <- prefLL + sufLL
  mode <- which.max(totalLL)
  
  # To also return fitted probabilities p_hat, we must actually compute the two
  # segment fits at the chosen mode (one PAVA on prefix, one on suffix).
  # (This is still O(n) total dominated by one extra pass.)
  pava_fit_binom <- function(yseg, wseg) {
    # same as blocks but return expanded fit
    n2 <- length(yseg)
    sw  <- numeric(n2); swy <- numeric(n2); lvl <- numeric(n2)
    st  <- integer(n2); en  <- integer(n2)
    m <- 0L
    for (i in seq_len(n2)) {
      m <- m + 1L
      sw[m]  <- wseg[i]
      swy[m] <- wseg[i] * yseg[i]
      lvl[m] <- swy[m] / sw[m]
      st[m] <- i; en[m] <- i
      while (m >= 2L && lvl[m-1L] > lvl[m]) {
        sw[m-1L]  <- sw[m-1L]  + sw[m]
        swy[m-1L] <- swy[m-1L] + swy[m]
        lvl[m-1L] <- swy[m-1L] / sw[m-1L]
        en[m-1L]  <- en[m]
        m <- m - 1L
      }
    }
    fit <- numeric(n2)
    for (j in seq_len(m)) fit[st[j]:en[j]] <- lvl[j]
    fit
  }
  
  left_fit  <- pava_fit_binom(y[1:mode], w[1:mode])
  right_fit <- rev(pava_fit_binom(rev(y[mode:n]), rev(w[mode:n])))  # nonincreasing on mode:n
  
  p_hat <- numeric(n)
  p_hat[1:mode] <- left_fit
  p_hat[mode:n] <- right_fit  # overwrites p_hat[mode] consistently for the mode
  
  list(
    p = p_hat,
    mode_index = mode,
    mode_L = L[mode],
    logLik = totalLL[mode],
    prefLL = prefLL,
    sufLL = sufLL
  )
}


# 3 Function for generating data according to ILC-IRT ----

# Function for generating data for the kth condition
gen_data <- function(k, params, seed){
  
  # 1) SETUP
  
  # Filter out other conditions than the kth condition from the parameter data frame
  fun_par <- params[params$sim == k,]
  
  # Extract all parameters
  N <- fun_par$N
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
  
  # 2) GENERATE THE SAMPLE
  
  # Set seed
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # Generate person parameters
  sample <- data.frame(mvrnorm(n = N, mu = mu_P, Sigma = Sigma_P, empirical = F))
  colnames(sample) <- c("phi", "theta", "tau")
  
  # Generate person-item engagement, correctness, and log RTs
  sample$d <- numeric(nrow(sample))
  sample$u <- numeric(nrow(sample))
  sample$l <- numeric(nrow(sample))
  for(i in 1:N){
    
    # Engagement indicator
    sample$d[i] <- ifelse(p_engagement(phi = sample$phi[i], iota = iota) < runif(1), 0, 1)
    
    # Response correctness
    if(sample$d[i] == 1){
      sample$u[i] <- ifelse(p_correct(theta = sample$theta[i], 
                                      a = a, 
                                      b = b, 
                                      c = c)
                            < runif(1), 0, 1)
    }else{
      sample$u[i] <- runif(1) < g
    }
    
    # Log RT
    if(sample$d[i] == 1){
      sample$l[i] <- rnorm(1, mean = beta - sample$tau[i], sd = alpha)
    }else{
      sample$l[i] <- rnorm(1, mean = mu_c, sd = sigma_c)
    }
  }
  
  # Reorder data w.r.t. ascending log RT
  ind <- order(sample$l)
  sample <- sample[ind,]
  
  # Done
  sample
}

# 4 Functions for computing the posterior probability of being in the engaged state given the observed data ----

# We utilised OpenAI's GPT-5.2 model to write the following functions

# Function for creating a "pre-computation object" for model parameters
# This computes everything that does not depend on l0
make_precomp <- function(n_gh = 25, y_cor = T, iota, a, b, c, alpha, beta, mu_P, Sigma_P){
  
  # precompute linear algebra once
  S11 <- Sigma_P[1,1]
  sd_phi <- sqrt(S11)
  A21 <- Sigma_P[2,1] / S11
  sd_th_cphi <- sqrt(Sigma_P[2,2] - Sigma_P[2,1]^2 / S11)
  
  S12 <- Sigma_P[1:2, 1:2]
  S12_inv <- solve(S12)
  B <- Sigma_P[3, 1:2] %*% S12_inv              # 1x2
  sd_tau_c <- sqrt(Sigma_P[3,3] - Sigma_P[3,1:2] %*% S12_inv %*% Sigma_P[1:2,3])
  sd_L <- sqrt(alpha^2 + sd_tau_c^2)

  # Gauss–Hermite:
  i <- 1:(n_gh-1)
  J <- matrix(0, n_gh, n_gh)
  J[cbind(i, i+1)] <- sqrt(i/2)
  J[cbind(i+1, i)] <- sqrt(i/2)
  eig <- eigen(J, symmetric = TRUE)
  x <- eig$values
  w <- (eig$vectors[1, ]^2) * sqrt(pi)
  o <- order(x)
  x <- x[o]
  w <- w[o]
  
  # build 2D grid
  X1 <- rep(x, times = n_gh)
  X2 <- rep(x, each  = n_gh)
  W  <- rep(w, times = n_gh) * rep(w, each = n_gh)
  
  # Convert to Z ~ N(0,1): z = sqrt(2) * x
  z1 <- sqrt(2) * X1
  z2 <- sqrt(2) * X2
  
  # transform to (phi, theta)
  phi <- mu_P[1] + sd_phi * z1
  mu_theta <- mu_P[2] + A21 * (phi - mu_P[1])
  theta <- mu_theta + sd_th_cphi * z2
  
  # mu_tau(phi,theta) = mu3 + B %*% ([phi,theta]-mu12)
  mu_tau <- as.numeric(mu_P[3] + B[1,1]*(phi - mu_P[1]) + B[1,2]*(theta - mu_P[2]))
  
  # h(phi,theta) part (no l0)
  if(y_cor){
    h <- p_engagement(phi, iota = iota) * p_correct(theta, a = a, b = b, c = c)
  }else{
    h <- p_engagement(phi, iota = iota) * (1-p_correct(theta, a = a, b = b, c = c))
  }
  
  # Done
  list(
    W = W,
    h = h,
    mu_tau = mu_tau,
    sd_L = sd_L,
    beta = beta
  )
}

# Function for computing the integral over a vector of log RTs (l0) 
# using the pre-computation object
integrate_over_l0 <- function(l0, pre){
  
  # compute density matrix: M x K via outer
  dens <- outer(pre$beta - l0, pre$mu_tau, function(xx, mm)
    dnorm(xx, mean = mm, sd = pre$sd_L)
  )
  
  as.numeric((dens %*% (pre$W * pre$h))/pi)
}

# Function for computing the posterior P(D = 1 | l, y) for each observation
compute_p_d1_posterior <- function(k, params, sim_data, n_gh = 25, p_d1_df){
  
  # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #
  # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #
  
  # 1) SETUP
  
  # Filter out other conditions from the parameter data frame
  fun_par <- params[params$sim == k,]
  
  # Extract all parameters
  N <- fun_par$N
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
  
  # Extract pre-computed P(Delta = 1 | iota) for this condition
  p_d1 <- p_d1_df$p_d1[p_d1_df$iota == iota]
  
  # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #
  # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #
  
  # 2) COMPUTE POSTERIOR FOR MR FOR EACH OBSERVED Y AND L
  
  # Create precomputation objects for the parameters
  
  # Correct answers
  pre_u1 <- make_precomp(
    n_gh = n_gh,
    y_cor = T,
    iota = iota, a = a, b = b, c = c,
    alpha = alpha, beta = beta,
    mu_P = mu_P, Sigma_P = Sigma_P)
  
  # Incorrect answers
  pre_u0 <- make_precomp(
    n_gh = n_gh,
    y_cor = F,
    iota = iota, a = a, b = b, c = c,
    alpha = alpha, beta = beta,
    mu_P = mu_P, Sigma_P = Sigma_P)
  
  # P(L = L, Y = y, Delta = d) for all y and d
  p_lu1d1 <- integrate_over_l0(sim_data$l, pre_u1)
  p_lu0d1 <- integrate_over_l0(sim_data$l, pre_u0)
  p_lu1d0 <- (1-p_d1)*dnorm(sim_data$l, mu_c, sigma_c)*g
  p_lu0d0 <- (1-p_d1)*dnorm(sim_data$l, mu_c, sigma_c)*(1-g)
  
  # P(L = l, Y = y, D = 1) and P(L = l, Y = y, D = 0)
  p_lyd1 <- sim_data$u*p_lu1d1 + (1-sim_data$u)*p_lu0d1
  p_lyd0 <- sim_data$u*p_lu1d0 + (1-sim_data$u)*p_lu0d0
  
  # P(L = l, Y = y)
  p_ly <- p_lyd1 + p_lyd0
  
  # P(D=1|L=l,Y=y)
  p_d1_cond_ly <- p_lyd1 / p_ly
  
  # done
  p_d1_cond_ly
}

# 5 Functions for computing thresholds ----

# Function for computing a MaxMI-22 threshold at the population level 
thr_pop_lvl_MaxMI22_f <- function(sim, params){
  
  # Filter out other items from the parameter data frame
  fun_par <- params[params$sim == sim,]
  
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
  
  # Done
  thr
}

# Function for computing a threshold by minimising the posterior P(Delta = 1 | l, y)
thr_real_lvl_MinMR_f <- function(sim_data_k, p_d1_cond_ly_k){
  
  # Extract log RT and engagement indicator
  l <- sim_data_k$l
  d <- sim_data_k$d
  
  # Sample size
  n <- length(l)
  
  # Compute MR based on sample
  mr_s <- cumsum(d) / n + (sum(1-d) - cumsum(1-d))/ n
  
  # Compute MR based on P(Delta = 1 | l, y)
  mr_p <- cumsum(p_d1_cond_ly_k) + c(rev(cumsum(rev(1-p_d1_cond_ly_k[-1]))), 0)
  
  # Compute a threshold by minimising MR based on P(Delta = 1 | l, y)
  thr <- l[dplyr::near(mr_p, min(mr_p))]
  
  # Also, compute MR at the threshold
  mr_s_at_thr <- mr_s[dplyr::near(mr_p, min(mr_p))] # based on sample
  mr_p_at_thr <- mr_p[dplyr::near(mr_p, min(mr_p))]/n # based on posterior, also scale it to be between 0 and 1
  
  # Done
  data.frame(thr = thr, mr_s = mr_s_at_thr, mr_p = mr_p_at_thr)
}

# Function for computing a MaxMI-22 threshold
thr_real_lvl_MaxMI22_f <- function(sim_data_k, p_d1_cond_ly_k, n_cut = 0, limit = T){
  
  # 1) Extract engagement indicator, log RT and response vectors
  l <- sim_data_k$l
  u <- sim_data_k$u
  d <- sim_data_k$d
  
  # 2) Compute the empirical conditional probabilities r = P(Y = 1 | X = 1) and 
  # q = P(Y = 1 | X = 0) and the empirical MIRT
  n <- length(l) # sample size
  p_x0 <- 1 - (1:n)/n # P(X=0)=1-P(X=1), needed for q = P(Y = 1 | X = 0)
  p_x0y1 <- mean(u)-cumsum(u)/n # P(X=0,Y=1)=P(Y=1)-P(Y=1,X=1), needed for q = P(Y = 1 | X = 0)
  r <- emp_cump_f(u) # r = P(Y=1|X=1), i.e., the CUMP curve
  q <- p_x0y1 / pmax(p_x0, .Machine$double.xmin) # q = P(Y=1|X=0)=P(X=0,Y=1)/P(X=0)
  mirt <- emp_mirt_x2y2_f(u)
  
  # 3) Compute misclassification rates based on sample and the posterior P(Delta = 1 | l,y)
  
  # Compute MR based on sample
  mr_s <- cumsum(d) / n + (sum(1-d) - cumsum(1-d))/ n
  
  # Compute MR based on P(Delta = 1 | l, y)
  mr_p <- cumsum(p_d1_cond_ly_k) + c(rev(cumsum(rev(1-p_d1_cond_ly_k[-1]))), 0)
  
  # 4) Cut first n_cut observations from the relevant vectors
  r0 <- r[(n_cut+1):n]
  q0 <- q[(n_cut+1):n]
  l0 <- l[(n_cut+1):n]
  mirt0 <- mirt[(n_cut+1):n]
  mr_s0 <- mr_s[(n_cut+1):n]
  mr_p0 <- mr_p[(n_cut+1):n]
  
  # 5) If limit = T, locate the first intersection of r and q, i.e., the first time point 
  # where rq = r - q switches its sign
  if(limit){
    rq0 <- r0 - q0
    l_rq <- NA
    for(m in 1:(length(rq0)-1)){ # skip the last observation
      if(rq0[m] == 0){
        l_rq <- l0[m]
        break()
      }else if(rq0[m] < 0 & rq0[m+1] > 0 | rq0[m] > 0 & rq0[m+1] < 0){
        l_rq <- (l0[m] + l0[m+1])/2
        break()
      }
    }
  }else{
    l_rq <- NA
  }
  
  # 5) Compute the MaxMI-22 threshold
  
  # If l_rq was found, find the global maximum of empirical MIRT between the minimum log RT in the sample and l_rq
  if(!is.na(l_rq)){
    
    # Limit the search space:
    l1 <- l0[l0 < l_rq]
    mirt1 <- mirt0[l0 < l_rq]
    mr_s1 <- mr_s0[l0 < l_rq]
    mr_p1 <- mr_p0[l0 < l_rq]
    
    # Compute threshold:
    thr <- l1[which.max(mirt1)]
    
    # Also, compute MR at the threshold
    mr_s_at_thr <- mr_s1[which.max(mirt1)] # based on sample
    mr_p_at_thr <- mr_p1[which.max(mirt1)]/n # based on posterior, also scale it between 0 and 1
    
  # If l_rq was not found, find the global maximum in the whole sample
  }else{
    # Compute threshold:
    thr <- l0[which.max(mirt0)]
    
    # Also, compute MR at the threshold
    mr_s_at_thr <- mr_s0[which.max(mirt0)] # based on sample
    mr_p_at_thr <- mr_p0[which.max(mirt0)]/n # based on posterior, also scale it between 0 and 1
  }
  
  # Done
  res <- list()
  res$l_rq <- l_rq
  res$thr_df <- data.frame(thr = thr, mr_s = mr_s_at_thr, mr_p = mr_p_at_thr)
  res
}

# Function for computing a MaxMI-22 threshold where MIRT is first smoothed with unimodal monotone regression
thr_real_lvl_MaxMI22_smooth_f <- function(sim_data_k, p_d1_cond_ly_k, n_cut = 0, limit = T, eps = .Machine$double.xmin){
  
  # 1) Extract engagement indicator, log RT and response vectors
  l <- sim_data_k$l
  u <- sim_data_k$u
  d <- sim_data_k$d
  
  # sample size
  n <- length(l)
  
  # 2) Smoothen the empirical MIRT based on unimodal monotone regression
  
  # 2.1) Fit unimodal monotone regression with mode and antimode on the response vector u
  p1 <- unimodal_regression(l, u)$p # one mode
  p2 <- 1 - unimodal_regression(l, 1 - u)$p # one antimode
  
  # 2.2) Compute log-likelihoods of the fits
  loglik1 <- sum(u * log(pmax(p1, eps)) + (1 - u) * log(pmax(1 - p1, eps)))
  loglik2 <- sum(u * log(pmax(p2, eps)) + (1 - u) * log(pmax(1 - p2, eps)))
  
  # 2.3) Compute smoothed r = P(U = 1 | X = 1), q = P(U = 1 | X = 0), and MIRT
  
  # Choose the fit with one mode if loglik1 >= loglik2
  if(loglik1 >= loglik2){
    p_x1 <- (1:n)/n
    r <- (mean(u)/(p_x1+eps))*cumsum(p1)/sum(p1)
    q <- (mean(u)/(1-p_x1+eps))*(1-cumsum(p1)/sum(p1))
    mirt <- h_f(mean(u)) - (1-p_x1)*h_f(q) - p_x1*h_f(r)
    
  # Choose fit with one antimode if loglik1 < loglik2
  }else{ 
    p_x1 <- (1:n)/n
    r <- (mean(u)/(p_x1+eps))*cumsum(p2)/sum(p2)
    q <- (mean(u)/(1-p_x1+eps))*(1-cumsum(p2)/sum(p2))
    mirt <- h_f(mean(u)) - (1-p_x1)*h_f(q) - p_x1*h_f(r)
  }
  
  # 3) Compute misclassification rates based on sample and the posterior P(Delta = 1 | l,y)
  
  # Compute MR based on sample
  mr_s <- cumsum(d) / n + (sum(1-d) - cumsum(1-d))/ n
  
  # Compute MR based on P(Delta = 1 | l, y)
  mr_p <- cumsum(p_d1_cond_ly_k) + c(rev(cumsum(rev(1-p_d1_cond_ly_k[-1]))), 0)
  
  # 4) Cut first n_cut observations from the relevant vectors
  r0 <- r[(n_cut+1):n]
  q0 <- q[(n_cut+1):n]
  l0 <- l[(n_cut+1):n]
  mirt0 <- mirt[(n_cut+1):n]
  mr_s0 <- mr_s[(n_cut+1):n]
  mr_p0 <- mr_p[(n_cut+1):n]
  
  # 5) If limit = T, locate the first intersection of r and q, i.e., the first time point 
  # where rq = r - q switches its sign
  if(limit){
    rq0 <- r0 - q0
    l_rq <- NA
    for(m in 1:(length(rq0)-1)){ # skip the last observation
      if(rq0[m] == 0){
        l_rq <- l0[m]
        break()
      }else if(rq0[m] < 0 & rq0[m+1] > 0 | rq0[m] > 0 & rq0[m+1] < 0){
        l_rq <- (l0[m] + l0[m+1])/2
        break()
      }
    }
  }else{
    l_rq <- NA
  }
  
  # 5) Compute the MaxMI-22 threshold
  
  # If l_rq was found, find the global maximum of empirical MIRT between the minimum log RT in the sample and l_rq
  if(!is.na(l_rq)){
    
    # Limit the search space:
    l1 <- l0[l0 < l_rq]
    mirt1 <- mirt0[l0 < l_rq]
    mr_s1 <- mr_s0[l0 < l_rq]
    mr_p1 <- mr_p0[l0 < l_rq]
    
    # Compute threshold:
    thr <- l1[which.max(mirt1)]
    
    # Also, compute MR at the threshold
    mr_s_at_thr <- mr_s1[which.max(mirt1)] # based on sample
    mr_p_at_thr <- mr_p1[which.max(mirt1)]/n # based on posterior, also scale it between 0 and 1
    
    # If l_rq was not found, find the global maximum in the whole sample
  }else{
    # Compute threshold:
    thr <- l0[which.max(mirt0)]
    
    # Also, compute MR at the threshold
    mr_s_at_thr <- mr_s0[which.max(mirt0)] # based on sample
    mr_p_at_thr <- mr_p0[which.max(mirt0)]/n # based on posterior, also scale it between 0 and 1
  }
  
  # Done
  res <- list()
  res$umr_loglik <- list("mode" = loglik1, "antimode" = loglik2)
  res$l_rq <- l_rq
  res$thr_df <- data.frame(thr = thr, mr_s = mr_s_at_thr, mr_p = mr_p_at_thr)
  res
}

# Function for computing a CUMP threshold
thr_real_lvl_CUMP_f <- function(sim_data_k, p_d1_cond_ly_k, g=.25, n_cut = 0){
  
  # 1) Extract engagement indicator, log RT and response vectors
  l <- sim_data_k$l
  u <- sim_data_k$u
  d <- sim_data_k$d
  
  # sample size
  n <- length(l)
  
  # 2) Compute CUMP
  cump <- emp_cump_f(u)
  
  # 3) Compute misclassification rates based on sample and the posterior P(Delta = 1 | l,y)
  
  # Compute MR based on sample
  mr_s <- cumsum(d) / n + (sum(1-d) - cumsum(1-d))/ n
  
  # Compute MR based on P(Delta = 1 | l, y)
  mr_p <- cumsum(p_d1_cond_ly_k) + c(rev(cumsum(rev(1-p_d1_cond_ly_k[-1]))), 0)
  
  # 4) Cut first n_cut observations from the relevant vectors
  cump0 <- cump[(n_cut+1):n]
  l0 <- l[(n_cut+1):n]
  mr_s0 <- mr_s[(n_cut+1):n]
  mr_p0 <- mr_p[(n_cut+1):n]
  
  # 5) Compute the threshold
  thr <- ifelse(mean(u) > g,
                tail(l0[cump0 <= g], n = 1), 
                tail(l0[cump0 >= g], n = 1))
  
  # Also, compute MR at the threshold
  mr_s_at_thr <- tail(mr_s0[l0 <= thr], 1) # based on sample
  mr_p_at_thr <- tail(mr_p0[l0 <= thr], 1)/n # based on posterior, also scale it between 0 and 1
  
  # Done
  data.frame(thr = thr, mr_s = mr_s_at_thr, mr_p = mr_p_at_thr)
}

# Function for computing an MLN threshold
thr_real_lvl_MLN_f <- function(sim_data_k, p_d1_cond_ly_k, seed_k){
  
  # 1) Extract engagement indicator, log RT and response vectors
  l <- sim_data_k$l
  u <- sim_data_k$u
  d <- sim_data_k$d
  
  # sample size
  n <- length(l)
  
  # 2) Compute misclassification rates based on sample and the posterior P(Delta = 1 | l,y)
  
  # Compute MR based on sample
  mr_s <- cumsum(d) / n + (sum(1-d) - cumsum(1-d))/ n
  
  # Compute MR based on P(Delta = 1 | l, y)
  mr_p <- cumsum(p_d1_cond_ly_k) + c(rev(cumsum(rev(1-p_d1_cond_ly_k[-1]))), 0)
  
  # 3) Fit a two component model to the data assuming normal distributions for log RT
  set.seed(seed_k)
  mixmod2_error <- FALSE
  tryCatch(mixmod2 <- normalmixEM(x = l, k = 2),
           error = function(e) {mixmod2_error <<- TRUE})
  
  # 4) Determine whether the estimation converged.
  
  # If an error occured, the mixmod2 won't exist in the environment for this function
  if(mixmod2_error){
    mixmod2_no_convergence <- NA
    
  # If an error did not occur, we can determine convergence from the number of iterations
  # (Max iteration number is 1000)
  }else{
    mixmod2_no_convergence <- length(mixmod2$all.loglik) == 1001
  }
  
  # 5) Compute BIC values for one- and two-component models
  
  # If an error occurred when fitting the two-component model, 
  # or the model estimation did not converge, 
  # we do not compute BIC for the two-component model
  if(mixmod2_error | (!is.na(mixmod2_no_convergence) & mixmod2_no_convergence)){
    
    # Log-likelihoods
    loglik1 <- (-n/2)*(log(2*pi*((n-1)/n)*sd(l)^2)+1) # one-component model's log-likelihood
    loglik2 <- NA # no log-likelihood for the two-component model in this case
    
    # Calculate BIC for both models
    BIC1 <- My_BIC(loglik = loglik1, q = 2, n = n)
    BIC2 <- NA # No likelihood for the two-component model, so BIC2 = NA
    
    # Save the two-component mixture model parameters as missing values
    mixmod_params <- list("lambda" = NA, "mu" = NA, "sigma" = NA)
    
  # If the estimation converged without error, we compute BIC for the two-component model
  }else{
    
    # Log-likelihoods
    loglik1 <- (-n/2)*(log(2*pi*((n-1)/n)*sd(l)^2)+1) # one-component model's log-likelihood
    loglik2 <- mixmod2$loglik # # two-component model's log-likelihood
    
    # Calculate BIC for both models
    BIC1 <- My_BIC(loglik = loglik1, q = 2, n = n)
    BIC2 <- My_BIC(loglik = loglik2, q = 5, n = n)
    
    # Save the MMA objects
    mixmod_params <- list("lambda" = mixmod2$lambda,
                          "mu" = mixmod2$mu, 
                          "sigma" = mixmod2$sigma)
  }
  
  # Save the computed objects
  number_of_components_analysis <- list("LL1" = loglik1, "LL2" = loglik2,
                                        "BIC1" = BIC1, "BIC2" = BIC2, 
                                        "mixmod2_no_error" = !mixmod2_error, 
                                        "mixmod2_converged" = !mixmod2_no_convergence)
  
  # 6) Determine whether to continue with threshold estimation or not
  
  # The conditions are...
  # a) The two-component model must fit the data better than the one-component model
  # b) The parameter estimates must correspond to our conceptual understanding of the RT distribution when RRB is present
  ### (in practice, the RRB weight must be below 0.5 and the RT mode in the RRB component must be smaller than the RT mode in the SB component)
  # c) The interval between the RT modes must contain only one local minimum of the estimated RT density
  
  # Condition a)
  if(!is.na(BIC2) & BIC2 < BIC1){
    condA <- TRUE
  }else{
    condA <- FALSE
  }
  
  # Conditions b) and c)
  # As long as the model converged without error, we can check conditions b) and c)
  if(!is.na(BIC2)){
    
    # Extract mixture model parameter estimates
    
    # lambda, mu, sigma
    lambda_RRB <- mixmod2$lambda[which.min(mixmod2$mu)]
    mu_RRB <- mixmod2$mu[which.min(mixmod2$mu)]
    sigma_RRB <- mixmod2$sigma[which.min(mixmod2$mu)]
    mu_SB <- mixmod2$mu[which.max(mixmod2$mu)]
    sigma_SB <- mixmod2$sigma[which.max(mixmod2$mu)]
    
    # RT modes
    mode_RT_RRB <- mu_RRB - sigma_RRB^2
    mode_RT_SB <- mu_SB - sigma_SB^2
    
    # Condition b):
    if(lambda_RRB < 0.5 & mode_RT_RRB < mode_RT_SB){
      condB <- TRUE
    }else{
      condB <- FALSE
    }
    
    # Estimate RT density as a function of log RT
    p_t_est_f <- function(l) lambda_RRB*dnorm(l, mu_RRB, sigma_RRB)/exp(l)+(1-lambda_RRB)*dnorm(l, mu_SB, sigma_SB)/exp(l)
  
    # Determine roughly whether the RT density has a global minimum between the modes: 
    # We find all local minima between the RT modes in .0001 wide segments
    if(mode_RT_RRB < mode_RT_SB){
      l_seg <- seq(mode_RT_RRB, mode_RT_SB, .0001)
      loc_min_p_t <- local_minima(p_t_est_f(l_seg))
      
      # Condition c)
      if(length(loc_min_p_t) == 1){
        condC <- TRUE
      }else{
        condC <- FALSE
      }
    }else{
      condC <- FALSE
    }
    
  # If the model did not converge, conditions B and C are useless
  }else{
    condB <- FALSE
    condC <- FALSE
  }
  
  # 7) If all conditions held, we can proceed with threshold computation
  if(condA & condB & condC){
    
    # Compute a threshold
    thr <- optimize(p_t_est_f, interval = c(mode_RT_RRB, mode_RT_SB))$minimum
    
    # Also, compute MR at the threshold
    # If the threshold is smaller than the smallest log RT in the sample, set MR to the RRB proportion in sample, 
    if(thr < min(l)){
      mr_s_at_thr <- 1-mean(d) # based on sample
      mr_p_at_thr <- 1-mean(d)
      
    # If the threshold is greater than or equal to the smallest log RT in the sample, compute MR normally
    }else{
      mr_s_at_thr <- tail(mr_s[l <= thr], 1) # based on sample
      mr_p_at_thr <- tail(mr_p[l <= thr], 1)/n # based on posterior, also scale it between 0 and 1
    }
  }else{
    thr <- NA
    mr_s_at_thr <- NA
    mr_p_at_thr <- NA
  }
  
  # done
  res <- list()
  res$mixmod_params <- mixmod_params
  res$number_of_components_analysis <- number_of_components_analysis
  res$threshold_computation_conditions <- list(condA = condA, condB = condB, condC = condC)
  res$thr_df <- data.frame(thr = thr, mr_s = mr_s_at_thr, mr_p = mr_p_at_thr)
  res
}

# END OF SCRIPT ----