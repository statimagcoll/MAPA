# Function to compute one-step
# muhat: predicted y values
# y: actual y
# transform: estimator transformation (squared logit ("q") by default). Options: "none", "no" for no transform; "l" or "logit" for logit; "q" or "squared logit" for squared logit, "z", "f" or "fisher" for Fisher
#' @importFrom stats cor
#' @noRd
cor_os = function(muhat, y, transform = c("q")){
  p = cor(muhat, y)
  if (transform %in% c("none", "no")){
    r = p + mean(if_full(muhat, y))
  }
  if (transform %in% c("l", "logit")){
    p = ifelse(p < 0, 0, ifelse(p > 1, 1, p))
    r = expit(logit(p) + mean(if_logit(muhat, y)))
  }
  if (transform %in% c("q", "squared logit")){
    p = ifelse(p < -1, -1, ifelse(p > 1, 1, p))
    r = sqrt(expit(logit(p^2) + mean(if_logit_sq(muhat, y)))) # fixed from if_full to if_logit_sq 08/12/2025
  }
  if (transform %in% c("z", "f", "fisher")){
    p = ifelse(p < -1, -1, ifelse(p > 1, 1, p))
    r = tanh(atanh(p) + mean(if_fisher(muhat, y)))
  }

  return(r)
}

# Function to compute difference between two one-step estimators (os2 - os1) (with confidence interval)
# os2: one-step estimate to subtract from
# os1: one-step estimate to subtract
# influence2: vector of influence functions for parameter 2 (corresponding to os2) for entire sample
# influence1: vector of influence functions for parameter 1 (corresponding to os1) for entire sample
# transform: estimator transformation (squared logit ("q") by default). Options: "none", "no" for no transform; "l" or "logit" for logit; "q" or "squared logit" for squared logit, "z", "f" or "fisher" for Fisher
# cl: confidence level
#' @noRd
cor_os_diff = function(os2, os1, influence2, influence1, transform = "q", cl = 0.95){
  if(length(influence1) != length(influence2)){
    stop("Influence function estimates must come from the same sample")
  }
  n = length(influence1)

  # point estimate: difference
  point_est = os2 - os1

  # variance: v(os1) + v(os2) - 2*cov(os1, os2)
  # computed as sum of the empirical variances of the influence functions minus 2 times the empirical covariance of the influence function estimates
  v1 = var(influence1)
  v2 = var(influence2)
  cov12 = cov(influence1, influence2)
  if (transform %in% c("none", "no")){
    var_factor = function(rho_hat) 1
  }
  if (transform %in% c("l", "logit")){
    var_factor = function(rho_hat){
      rho_hat*(1-rho_hat)
    }
  }
  if (transform %in% c("q", "squared logit")){
    var_factor = function(rho_hat){
      0.5*rho_hat*(1-rho_hat^2)
    }
  }
  if (transform %in% c("z", "f", "fisher")){
    var_factor = function(rho_hat){
      (1-rho_hat^2)
    }
  }

  v_diff = v1*var_factor(os1)^2 + v2*var_factor(os2)^2 - 2*cov12*var_factor(os1)*var_factor(os2)

  lower = point_est - qnorm(cl+(1-cl)/2)*sqrt(v_diff/n)
  upper = point_est + qnorm(cl+(1-cl)/2)*sqrt(v_diff/n)

  return(c("estimate" = point_est, "lCI" = lower, "uCI" = upper))
}

# Function to compute ratio of Cohen's f^2 via two one-step (squared logit) estimators (with confidence interval)
# os2: one-step estimate (numerator)
# os1: one-step estimate (denominator)
# influence2: vector of influence functions for parameter 2 (corresponding to os2) for entire sample
# influence1: vector of influence functions for parameter 1 (corresponding to os1) for entire sample
# cl: confidence level
#' @noRd
cor_os_ratio = function(os2, os1, influence2, influence1, cl = 0.95){
  #browser()
  if(length(influence1) != length(influence2)){
    stop("Influence function estimates must come from the same sample")
  }
  n = length(influence1)

  # point estimate: exponeniated log ratio
  point_est = exp(logit(os2^2)-logit(os1^2))

  # variance: v(os1) + v(os2) - 2*cov(os1, os2)
  # computed as sum of the empirical variances of the influence functions minus 2 times the empirical covariance of the influence function estimates
  v1 = var(influence1)
  v2 = var(influence2)
  cov12 = cov(influence1, influence2)

  v_diff = v1 + v2 - 2*cov12

  lower = exp(log(point_est) - qnorm(cl+(1-cl)/2)*sqrt(v_diff/n))
  upper = exp(log(point_est) + qnorm(cl+(1-cl)/2)*sqrt(v_diff/n))

  return(c("estimate" = point_est, "lCI" = lower, "uCI" = upper))
}

# logit and expit functions
#' @noRd
logit = function(x){
  log(x/(1-x))
}

#' @noRd
expit = function(x){
  1/(1+exp(-x))
}

## Influence Functions ##
# Functions for empirical influence function (returns vector)

# Fully derived from either Pearson or Simplified parameter
#' @noRd
if_full = function(muhat, y){
  # components
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  # influence function
  return((2*y*(muhat - my) - muhat^2 + my^2)/(2*sqrt(vmuhat*vy)) -
           (sqrt(vmuhat)*(y - my)^2/(2*(vy)^(3/2))))
}



# IF(logit(rho_0))
#' @noRd
if_logit = function(muhat, y){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  # fixed 3/5/25
  return((sqrt(vy)/sqrt(vmuhat) + sqrt(vy)/(sqrt(vy) - sqrt(vmuhat))) * if_full(muhat, y))
}

# IF(logit(rho_0^2))
#' @noRd
if_logit_sq = function(muhat, y){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)

  return(((y-my)^2*(vy-vmuhat) - vy*(y-muhat)^2)/(vmuhat*(vy-vmuhat)))
}

#' @noRd
if_fisher = function(muhat, y){
  return((1 / (1 - cor(muhat, y)^2)) * if_full(muhat, y))
}


### AUC Functions ###
#' @noRd
if_auc_direct <- function(muhat, y, pi_hat, F0, F1, f0, f1, ...){
  auc_plugin <- auc(roc(y, muhat, ...))
  if_direct <- (y == 1)*(F0(muhat) - auc_plugin)/pi_hat +
    (y == 0)*(1 - F1(muhat) - auc_plugin)/(1 - pi_hat)
  return(if_direct)
}

# Additional contribution from perturbing mu, arguments similar to above
#' @noRd
if_auc_mu <- function(muhat, y, f0, f1){
  if_mu <- (y - muhat) * (f1(muhat) - f0(muhat))
  return(if_mu)
}

# Simple AUC os function that computes fold-specific os
#' @noRd
auc_os_fold <- function(muhat, y, pi_hat, F0, F1, f0, f1, ...){
  # calculate AUC plug-in
  auc_plugin <- auc(roc(y, muhat, ...))
  # calculate direct IF values
  if_direct <- if_auc_direct(muhat, y, pi_hat, F0, F1, f0, f1, ...)
  # calculate additional IF
  if_mu <- if_auc_mu(muhat, y, f0, f1)

  # OS = plug-in + mean(if_direct + if_mu)
  os <- auc_plugin + mean(if_direct + if_mu)
  return(os)
}

# Function to compute difference between two one-step estimators (os2 - os1) (with confidence interval)
# os2: one-step estimate to subtract from
# os1: one-step estimate to subtract
# influence2: vector of influence functions for parameter 2 (corresponding to os2) for entire sample
# influence1: vector of influence functions for parameter 1 (corresponding to os1) for entire sample
# cl: confidence level
#' @noRd
auc_os_diff = function(os2, os1, influence2, influence1, cl = 0.95){
  if(length(influence1) != length(influence2)){
    stop("Influence function estimates must come from the same sample")
  }
  n = length(influence1)

  # point estimate: difference
  point_est = os2 - os1

  # variance: v(os1) + v(os2) - 2*cov(os1, os2)
  # computed as sum of the empirical variances of the influence functions minus 2 times the empirical covariance of the influence function estimates
  v1 = var(influence1)
  v2 = var(influence2)
  cov12 = cov(influence1, influence2)

  v_diff = v1 + v2 - 2*cov12

  lower = point_est - qnorm((1+cl)/2)*sqrt(v_diff/n)
  upper = point_est + qnorm((1+cl)/2)*sqrt(v_diff/n)

  return(c("estimate" = point_est, "lCI" = lower, "uCI" = upper, "se" = sqrt(v_diff/n)))
}
