# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for user-defined functions that are used in the 
# "1.2_pisa2022_results.R" script.

# Author:
# Santeri Holopainen
# Turku Research Institute for Learning Analytics (TRILA)
# University of Turku

# Date: 11.9.2026

# R version: 4.5.0

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# Function for calculating a fitted model's Bayesian Information Criterion (BIC)
My_BIC <- function(loglik, q, n){
  # loglik = log-likelihood of the fitted model
  # q = number of parameters in the model
  # n = sample size
  q*log(n)-2*loglik
}

# Function for finding local minima
local_minima <- function(x) {
  # x = the variable of interest
  s <- sign(diff(x))
  mask <- s != 0
  k <- which(diff(s[mask]) == 2)
  i <- which(s != 0)
  (i[k] + i[k+1] + 1) / 2
}

# Function for computing the entropy of a binary variable with probability p
h_f <- function(p, eps = .Machine$double.xmin){
  # p = probability of success
  # eps = very small value which is effectively treated as a zero in the log2 function
  -p*log2(pmax(p,eps))-(1-p)*log2(pmax(1-p,eps))
}

# Function for computing the empirical CUMP curve, 
# i.e., P(Y=1|X=1) = the empirical conditional probability of a correct response given X = 1 (T <= t_k)
emp_cump_f <- function(y){ 
  # y = response correctness variable, which is assumed to be sorted w.r.t. ascending response time
  cumsum(y) / 1:length(y)
}

# Function for computing empirical MIRT when both X and Y are binary
# This function returns a vector of the empirical MIRT values
emp_mirt_x2y2_f <- function(y, eps = .Machine$double.xmin){
  # y = response correctness variable, which is assumed to be sorted w.r.t. ascending response time
  # eps = very small value which is effectively treated as a zero in the log2 function
  n <- length(y) # sample size
  p_y1 <- mean(y) # P(Y = 1)
  p_x1 <- 1:n/n # P(X = 1) = P(T <= T_k)
  p_x1y1 <- (cumsum(y)/n) # P(X = 1, Y = 1)
  p_x0y1 <- p_y1-p_x1y1 # P(X = 0, Y = 1)
  p_y1_cond_x1 <- p_x1y1 / p_x1  # P(Y = 1 | X = 1) 
  p_y1_cond_x0 <- p_x0y1 / pmax((1-p_x1), eps) # P(Y = 1 | X = 0)
  h_f(p_y1)-(p_x1*h_f(p_y1_cond_x1) + (1-p_x1)*h_f(p_y1_cond_x0)) # empirical MIRT
}

# Function for computing empirical MIRT when... 
# - X is binary and
# - Y is the raw response variable and possibly has more than two unique values (q).
# This function returns a data frame of the empirical probabilities and MIRT.
# We utilised OpenAI's GPT-5.2 model for defining this function.
emp_mirt_x2yq_f <- function(t, y, eps = .Machine$double.xmin){
  # t = response time variable, which is assumed to be in ascending order
  # y = response correctness variable, which is assumed to be sorted w.r.t. ascending response time
  # eps = very small value which is effectively treated as a zero in the log2 function
  
  # Sample size
  n <- length(t)
  
  # Mutate y to factor
  y <- factor(y)
  
  # Extract levels of y (1, 2, ..., q)
  y_lev <- levels(y)
  
  # Compute empirical marginal probabilities of X and Y:
  p_x1 <- 1:n/n # P(X=1)
  p_x0 <- 1-p_x1 # P(X=0)
  p_y <- prop.table(table(y)) # P(Y=y), y = 1, 2, ..., q
  
  # Go through the levels of y. For each level, compute...
  # empirical conditional probabilities of X given Y,
  # empirical joint probabilities of X and Y, and
  # the MI term corresponding to the level of y.
  out_list <- list()
  for(q in 1:length(y_lev)){
    
    yq <- y_lev[q]
    Tq <- t[y == yq] # T | Y = y_q
    p_yq <- as.numeric(p_y[yq]) # P(Y = y_q)
    p_x1_cond_yq <- ecdf(Tq)(t) # P(X = 1 | Y = y_q)
    p_x1yq <- p_yq * p_x1_cond_yq # P(X = 1, Y = y_q)
    p_x0yq <- p_yq - p_x1yq # P(X = 0, Y = y_q)
    
    out_list[[q]] <- data.frame(
      t = t,
      y = yq,
      p_x1 = p_x1,
      p_x0 = p_x0,
      p_y = p_yq,
      p_x1y = p_x1yq,
      p_x0y = p_x0yq,
      p_x1_cond_y = p_x1_cond_yq,
      p_x0_cond_y = 1-p_x1_cond_yq,
      I_xy = 
        p_x1yq * (log2(pmax(pmin(p_x1yq, 1-eps), eps)) - log2(pmax(pmin(p_x1*p_yq, 1-eps), eps)))+
        p_x0yq * (log2(pmax(pmin(p_x0yq, 1-eps), eps)) - log2(pmax(pmin(p_x0*p_yq, 1-eps), eps)))
    )
  }
  
  # Bind the results together to a data frame
  out <- do.call(rbind, out_list)
  # This is now in "long" format w.r.t. the levels of y
  
  # Pivot the data frame wide w.r.t the levels
  out <- 
    out %>%
    pivot_wider(id_cols = c(t, p_x1, p_x0), names_from = y, values_from = p_y:I_xy, names_sep = "")
  
  # Compute the empirical MIRT as the sum of the terms corresponding to each level of y:
  out$mirt <- apply(out[,str_detect(colnames(out), "I_")], 1, sum)
  
  # Done
  out
}

# Function for computing empirical MIRT when...
# - X has three categories and
# - Y denotes response correctness.
# This function returns a data frame of the empirical probabilities and MIRT.
# We utilised OpenAI's GPT-5.2 model for defining this function.
emp_mirt_x3y2_f <- function(t, y, eps = .Machine$double.xmin){
  # t = response time variable, which is assumed to be in ascending order
  # y = response correctness variable, which is assumed to be sorted w.r.t. ascending response time
  # eps = very small value which is effectively treated as a zero in the log2 function
  
  # Extract unique response time values
  u <- sort(unique(t))
  
  # Find all pairs (t1,t2) with t1 < t2
  idx <- which(outer(u, u, `<`), arr.ind = TRUE)
  t1 <- u[idx[, 1]]
  t2 <- u[idx[, 2]]
  
  # Compute P(X = x), x = 0,1,2
  F1 <- ecdf(t)(t1) # P(T <= t1)
  F2 <- ecdf(t)(t2) # P(T <= t2)
  p_x2 <- F1 # P(X=2) = P(T <= t1)
  p_x1 <- F2 - F1 # P(X=1) = P(t1 < T <= t2)
  p_x0 <- 1 - F2 # P(X=0) = P(T > t2)
  
  # Compute P(Y = y), y = 0,1
  p_y0 <- mean(y == 0)
  p_y1 <- mean(y == 1)
  
  # Define functions for P(X = x | Y = y), x = 0, 1, 2; y = 0, 1
  Fy0_fun <- ecdf(t[y == 0])
  Fy1_fun <- ecdf(t[y == 1])
  
  # Compute P(X = x, Y = y), x = 0, 1, 2; y = 0, 1
  Fy0_t1 <- p_y0 * Fy0_fun(t1)
  Fy0_t2 <- p_y0 * Fy0_fun(t2)
  Fy1_t1 <- p_y1 * Fy1_fun(t1)
  Fy1_t2 <- p_y1 * Fy1_fun(t2)
  p_x2y0 <- Fy0_t1
  p_x1y0 <- Fy0_t2 - Fy0_t1
  p_x0y0 <- p_y0 - Fy0_t2
  p_x2y1 <- Fy1_t1
  p_x1y1 <- Fy1_t2 - Fy1_t1
  p_x0y1 <- p_y1 - Fy1_t2
  
  # Finally, compute empirical MIRT
  mirt <- 
    p_x0y0 * (log2(pmax(pmin(p_x0y0, 1-eps), eps)) - log2(pmax(pmin(p_x0*p_y0, 1-eps), eps)))+
    p_x1y0 * (log2(pmax(pmin(p_x1y0, 1-eps), eps)) - log2(pmax(pmin(p_x1*p_y0, 1-eps), eps)))+
    p_x2y0 * (log2(pmax(pmin(p_x2y0, 1-eps), eps)) - log2(pmax(pmin(p_x2*p_y0, 1-eps), eps)))+
    p_x0y1 * (log2(pmax(pmin(p_x0y1, 1-eps), eps)) - log2(pmax(pmin(p_x0*p_y1, 1-eps), eps)))+
    p_x1y1 * (log2(pmax(pmin(p_x1y1, 1-eps), eps)) - log2(pmax(pmin(p_x1*p_y1, 1-eps), eps)))+
    p_x2y1 * (log2(pmax(pmin(p_x2y1, 1-eps), eps)) - log2(pmax(pmin(p_x2*p_y1, 1-eps), eps)))
  
  # Done.
  data.frame(
    t1 = t1, 
    t2 = t2,
    p_y0 = p_y0,
    p_y1 = p_y1,
    p_x0 = p_x0,
    p_x1 = p_x1,
    p_x2 = p_x2,
    p_x0y0 = p_x0y0,
    p_x1y0 = p_x1y0,
    p_x2y0 = p_x2y0,
    p_x0y1 = p_x0y1,
    p_x1y1 = p_x1y1,
    p_x2y1 = p_x2y1,
    mirt = mirt
  )
}


# Function for computing a MaxMI-22 threshold
thr_MaxMI_x2y2_f <- function(d_k, n_cut = 0, limit = T){
  # d_k = data for item k
  # n_cut = number of observations with the smallest response times to cut from the data prior to threshold computation
  # In some cases, it is beneficial to clean the data a little bit from spurious noise
  # limit = should the search space for the maximum of the empirical MIRT be limited?
  
  # 1) Extract response time (RT) and response variables 
  # and sort them w.r.t ascending RT
  i <- order(d_k$Time)
  t <- d_k$Time[i]
  y <- d_k$Correctness[i]
  
  # 2) Compute the empirical conditional probabilities r = P(Y = 1 | X = 1) and 
  # q = P(Y = 1 | X = 0) and the empirical MIRT
  n <- length(t) # sample size
  p_x0 <- 1 - (1:n)/n # P(X=0)=1-P(X=1), needed for q = P(Y = 1 | X = 0)
  p_x0y1 <- mean(y)-cumsum(y)/n # P(X=0,Y=1)=P(Y=1)-P(Y=1,X=1), needed for q = P(Y = 1 | X = 0)
  r <- emp_cump_f(y) # r = P(Y=1|X=1), i.e., the CUMP curve
  q <- p_x0y1 / pmax(p_x0, .Machine$double.xmin) # q = P(Y=1|X=0)=P(X=0,Y=1)/P(X=0)
  mirt <- emp_mirt_x2y2_f(y)
  
  # 3) Cut first n_cut observations from the relevant vectors
  r0 <- r[(n_cut+1):n]
  q0 <- q[(n_cut+1):n]
  t0 <- t[(n_cut+1):n]
  mirt0 <- mirt[(n_cut+1):n]
  
  # 4) If limit = T, locate the first intersection of r and q, i.e., the first time point 
  # where rq = r - q switches its sign
  if(limit){
    rq0 <- r0 - q0
    t_rq <- NA
    for(m in 1:(length(rq0)-1)){ # skip the last observation
      if(rq0[m] == 0){
        t_rq <- t0[m]
        break()
      }else if(rq0[m] < 0 & rq0[m+1] > 0 | rq0[m] > 0 & rq0[m+1] < 0){
        t_rq <- (t0[m] + t0[m+1])/2
        break()
      }
    }
  }else{
    t_rq <- NA
  }
  
  # 5) Compute the MaxMI-22 threshold

  # If t_rq was found, find the global maximum of empirical MIRT between 0 and t_rq
  if(!is.na(t_rq)){
    t1 <- t0[t0 < t_rq]
    mirt1 <- mirt0[t0 < t_rq]
    thr <- t1[which.max(mirt1)]

  # If t_rq was not found, find the global maximum in the whole sample
  }else{
    thr <- t0[which.max(mirt0)]
  }
  
  # Done
  res <- list()
  res$t_rq <- t_rq
  res$thr <- thr
  res
}

# Function for computing a MaxMI-32 threshold
thr_MaxMI_x3y2_f <- function(d_k){
  # d_k = data for item k
  
  # 1) Extract response time (RT) and response variables 
  # and sort them w.r.t ascending RT
  i <- order(d_k$Time)
  t <- d_k$Time[i]
  y <- d_k$Correctness[i]
  
  # 2) Compute empirical MIRT when X has three categories and Y is binary
  # (recall that this results in a data frame where MIRT is in one column)
  d_x3y2 <- emp_mirt_x3y2_f(t = t, y = y)
  
  # 3) Compute the MaxMI-32 thresholds
  thr1 <- d_x3y2$t1[which.max(d_x3y2$mirt)]
  thr2 <- d_x3y2$t2[which.max(d_x3y2$mirt)]
  
  # Done
  res <- list()
  res$thr1 <- thr1
  res$thr2 <- thr2
  res
}

# Function for computing a MaxMI-2Q threshold
thr_MaxMI_x2yq_f <- function(d_k, n_cut = 0){
  # d_k = data for item k
  # n_cut = number of observations with the smallest response times to cut from the data prior to threshold computation
  # In some cases, it is beneficial to clean the data a little bit from spurious noise
  
  # 1) Extract response time (RT) and response variables 
  # and sort them w.r.t ascending RT
  i <- order(d_k$Time)
  t <- d_k$Time[i]
  y <- d_k$Response[i]
  
  # 2) Compute empirical MIRT when X is binary and Y is the raw response
  # (recall that this results in a data frame where MIRT is in one column)
  d_x2yq <- emp_mirt_x2yq_f(t = t, y = y)
  
  # 3) Cut first n_cut observations from the relevant vectors
  t0 <- d_x2yq$t[(n_cut+1):nrow(d_x2yq)]
  mirt0 <- d_x2yq$mirt[(n_cut+1):nrow(d_x2yq)]
  
  # 4) Compute the MaxMI-2Q threshold
  thr <- t0[which.max(mirt0)]
  
  # done
  thr
}

# Function for computing a CUMP threshold
thr_CUMP_f <- function(d_k, g_k, n_cut = 0){
  # d_k = data for item k
  # g_k = chance level for item k
  # n_cut = number of observations with the smallest response times to cut from the data prior to threshold computation
  ### In some cases, it is beneficial to clean the data a little bit from spurious noise
  
  # 1) Extract response time (RT) and response variables 
  # and sort them w.r.t ascending RT
  i <- order(d_k$Time)
  t <- d_k$Time[i]
  y <- d_k$Correctness[i]
  
  # 2) Compute CUMP
  cump <- emp_cump_f(y)
  
  # 3) Cut first n_cut observations from the relevant vectors
  t0 <- t[(n_cut+1):length(t)]
  cump0 <- cump[(n_cut+1):length(cump)]
  
  # 4) Compute the threshold
  thr <- ifelse(mean(y) > g_k,
                tail(t0[cump0 <= g_k], n = 1), 
                tail(t0[cump0 >= g_k], n = 1))
  
  # Done
  thr
}

# Function for computing an MLN threshold
thr_MLN_f <- function(d_k, seed_k){
  # d_k = data for item k
  # seed_k = seed number for reproducibility
  
  # 1) Extract response time (RT) and sort it in ascending order
  t <- sort(d_k$Time)
  l <- log(t) # log RT for fitting a mixture model
  n <- nrow(d_k) # sample size
  
  # 2) Fit a two component model to the data assuming normal distributions for log RT
  set.seed(seed_k)
  mixmod2_error <- FALSE
  tryCatch(mixmod2 <- normalmixEM(x = l, k = 2),
           error = function(e) {mixmod2_error <<- TRUE})
  
  # 3) Determine whether the estimation converged.
  
  # If an error occured, the mixmod2 won't exist in the environment for this function
  if(mixmod2_error){
    mixmod2_no_convergence <- NA
    
  # If an error did not occur, we can determine convergence from the number of iterations
  # (Max iteration number is 1000)
  }else{
    mixmod2_no_convergence <- length(mixmod2$all.loglik) == 1001
  }
  
  # 4) Compute BIC values for one- and two-component models
  
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
  
  # 5) Determine whether to continue with threshold estimation or not
  
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
  
  # 6) If all conditions held, we can proceed with threshold computation
  if(condA & condB & condC){
    # Take exponent because we are in the logarithm scale
    thr <- exp(optimize(p_t_est_f, interval = c(mode_RT_RRB, mode_RT_SB))$minimum)
    
  # If even one of the conditions did not hold, we cannot proceed
  }else{
    thr <- NA
  }
  
  # done
  res <- list()
  res$mixmod_params <- mixmod_params
  res$number_of_components_analysis <- number_of_components_analysis
  res$threshold_computation_conditions <- list(condA = condA, condB = condB, condC = condC)
  res$thr <- thr
  res
}

# END OF SCRIPT ----