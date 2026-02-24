#' This function computes the one-step estimator from Jones, Gadiyar, and Vandekar (2025).
#' Vector arguments are accepted. If different length arguments are passed they are dealt with in the usual way of R.
#' @param y A numeric vector of continuous values representing the real outcome values.
#' @param y_hat A numeric vector of continuous values representing the predicted outcome values.
#' @param folds A numeric vector representing the fold number of the corresponding y_hat and y value.
#' @param y_hat2 A numeric vector of continuous values representing the predicted outcome values from a different model.
#' @param level A number between 0 and 1 representing the confidence level.
#' @return Returns a list with a data frame of the one-step estimator, the logit-scale one-step estimator value, and its lower bound, upper bound, and standard error; a vector of the one-step estimate values by fold for the first model; the sample number; and the confidence level. If y_hat2 is passed as an argument, the one-step estimate for y_hat2 and its difference/ratio to the one-step estimate for y_hat1 are included with confidence intervals.
#'
#' @details The formula for calculating the one-step estimator is:
#' \eqn{\hat{\rho}_{OS} = \sqrt{\text{expit}\!\left(\text{logit}\!\left({\hat{\rho}_P}^2\right)\right)+\frac{1}{n}\sum_{i=1}^{n}\left\{\frac{(Y_i -\bar{Y})^2}{\hat{\text{Var}}(\hat{\mu)}}-\frac{(Y_i - \hat{\mu}(X_{i}))^2 \, \hat{\text{Var}}(Y)}{(\hat{\text{Var}}(Y) - \hat{\text{Var}}(\hat{\mu}(X))) \,\hat{\text{Var}}(\hat{\mu}(X))}\right\}}}
#'
#' @details The formula for calculating the one-step estimator confidence interval is:
#'
#' \eqn{\sqrt{\text{expit}\!\left(\text{logit}\!\left({\hat{\rho}_{OS}}^2\right)\right)\pm z_{1 - {\alpha/2}}\sqrt{\hat{\text{Var}}\!\left\{\frac{(Y - \bar{Y})^2}{\hat{\text{Var}}(\hat{\mu})}-\frac{(Y - \hat{\mu}(X))^2 \, \hat{\text{Var}}(Y)}{(\hat{\text{Var}}(Y) - \hat{\text{Var}}(\hat{\mu}(X))) \,\hat{\text{Var}}(\hat{\mu}(X))}\right\}/ n}}}
#'
#' @export

mapaFolds <- function(y, y_hat, folds, y_hat2 = NULL, level = 0.95){
  # Checking that variable lengths match
  if(length(y_hat) != length(y) | length(y_hat) != length(folds) | length(y) != length(folds)){
    stop("y, y_hat, and folds must be of the same length.")
  }

  # Creating storage vector for results and variables
  os_folds = c()
  inf_folds = c()
  n = length(y)

  # Looping through folds
  for(fold in sort(unique(folds))){
    # Calculating and storing OS estimator
    os_res = cor_os(y_hat[which(folds == fold)], y[which(folds == fold)])
    os_folds = c(os_folds, os_res)

    # Calculating and storing the influence vector
    inf_res = if_logit_sq(y_hat[which(folds == fold)], y[which(folds == fold)])
    inf_folds = c(inf_folds, inf_res)
  }

  # Calculating the mean OS estimate across folds
  os_estimate = mean(os_folds, na.rm = T)

  # Calculating the logit-scale estimate across folds
  os_estimate_logit = logit(os_estimate^2)

  # Calculating the standard error of the OS estimate
  stand_error = sqrt(var(inf_folds)/n)

  # Calculating the bounds of the OS estimate
  LB_est = sqrt(expit(os_estimate_logit - qnorm((1 + level) / 2)*stand_error))
  UB_est = sqrt(expit(os_estimate_logit - qnorm((1 - level) / 2)*stand_error))

  # Denoting the fold number that each OS fold estimate belongs to
  names(os_folds) <- sort(unique(folds))

  # Creating a results data frame for output
  estimate = data.frame(est = os_estimate, est_logit = os_estimate_logit, LB = LB_est, UB = UB_est, se = stand_error)


  if(!is.null(y_hat2)){
    # Checking that variable lengths match
    if(length(y_hat2) != length(y) | length(y_hat2) != length(folds) | length(y) != length(folds)){
      stop("y, y_hat2, and folds must be of the same length.")
    }

    # Creating storage vector for results and variables
    os_folds2 = c()
    inf_folds2 = c()
    n = length(y)

    # Looping through folds
    for(fold in sort(unique(folds))){
      # Calculating and storing OS estimator
      os_res = cor_os(y_hat2[which(folds == fold)], y[which(folds == fold)])
      os_folds2 = c(os_folds2, os_res)

      # Calculating the influence vector and storing
      inf_res = if_logit_sq(y_hat2[which(folds == fold)], y[which(folds == fold)])
      inf_folds2 = c(inf_folds2, inf_res)
    }

    # Calculating the mean OS estimate across folds
    os_estimate2 = mean(os_folds2, na.rm = T)

    # Calculating the logit-scale estimate across folds
    os_estimate2_logit = logit(os_estimate2^2)

    # Calculating standard error
    stand_error2 = sqrt(var(inf_folds2)/n)

    # Calculating the OS estimate bounds across folds
    LB_est2 = sqrt(expit(os_estimate2_logit - qnorm((1 + level) / 2)* stand_error2))
    UB_est2 = sqrt(expit(os_estimate2_logit - qnorm((1 - level) / 2)*stand_error2))

    # Denoting the fold number that each OS fold estimate belongs to
    names(os_folds2) <- sort(unique(folds))

    # Creating data frame output
    estimate2 = data.frame(est = os_estimate2, est_logit = os_estimate2_logit, LB = LB_est2, UB = UB_est2, se = stand_error2)

    # Calculating OS difference - first yhat vs. second yhat
    diff = cor_os_diff(os_estimate, os_estimate2, inf_folds, inf_folds2)
    os_diff = data.frame(est = diff[1], est_logit = NA, LB = diff[2], UB = diff[3], se = diff[4])

    # Calculating the OS ratio - first yhat vs. second yhat
    ratio = cor_os_ratio(os_estimate, os_estimate2, inf_folds, inf_folds2)
    os_ratio = data.frame(est = ratio[1], est_logit = NA, LB = ratio[2], UB = ratio[3], se = ratio[4])

    # Combining all data into one
    est_total = rbind(estimate, estimate2, os_diff, os_ratio)

    # Assigning row names
    est_total$metric <- c("yhat1", "yhat2", "Diff", "Ratio")
    rownames(est_total) <- NULL


    # Combining all outputs into one list
    out = list(est = est_total, fold.est = os_folds, fold.est2 = os_folds2, n = n, conf.level = level) # data frame version
  }
  else{
    # Combining all outputs into one list when there is no yhat2
    out = list(est = estimate, fold.est = os_folds, n = n, conf.level = level)

    # Naming the metric being measured in estimate
    out$est$metric = "yhat1"
  }

  # Returning the data frame
  return(out)
}


#' Perform multiple sample splits for one-step correlation estimator
#'
#' This function computes the one-step estimator from Jones et al. (2025). Multiple sample splits are recommended for
#' replicability and reproducibility. The final estimate and CI bounds can be taken as the median of the results from this function.
#'
#' @param y A vector of continuous values representing the real outcome values by fold.
#' @param y_hat A matrix of continuous values representing the predicted outcome values by fold. Each column corresponds to a different sample split.
#' @param folds A matrix representing the fold number of the corresponding y_hat and y value. Each column corresponds to a different sample split.
#' @param y_hat2 A matrix of continuous values representing the predicted outcome values by fold from a different model. Each column corresponds to a different sample split.
#' @param level A number between 0 and 1 representing the confidence level.
#' @return Returns a list with a data frame of the one-step estimator and its bounds per split, a matrix of the one-step estimate by fold and split, the sample number, and the confidence level. If y_hat2 is passed as an argument, the one-step estimate for y_hat2 and difference between the one-step estimate for y_hat1 and y_hat2 will be provided.
#'
#' @details The formula for calculating the one-step estimator is:
#' \eqn{\hat{\rho}_{OS} = \sqrt{\text{expit}\!\left(\text{logit}\!\left({\hat{\rho}_P}^2\right)\right)+\frac{1}{n}\sum_{i=1}^{n}\left\{\frac{(Y_i -\bar{Y})^2}{\hat{\text{Var}}(\hat{\mu)}}-\frac{(Y_i - \hat{\mu}(X_{i}))^2 \, \hat{\text{Var}}(Y)}{(\hat{\text{Var}}(Y) - \hat{\text{Var}}(\hat{\mu}(X))) \,\hat{\text{Var}}(\hat{\mu}(X))}\right\}}}
#'
#' @details The formula for calculating the one-step estimator confidence interval is:
#'
#' \eqn{\sqrt{\text{expit}\!\left(\text{logit}\!\left({\hat{\rho}_{OS}}^2\right)\right)\pm z_{1 - {\alpha/2}}\sqrt{\hat{\text{Var}}\!\left\{\frac{(Y - \bar{Y})^2}{\hat{\text{Var}}(\hat{\mu})}-\frac{(Y - \hat{\mu}(X))^2 \, \hat{\text{Var}}(Y)}{(\hat{\text{Var}}(Y) - \hat{\text{Var}}(\hat{\mu}(X))) \,\hat{\text{Var}}(\hat{\mu}(X))}\right\}/ n}}}
#'
#' @export
sampleSplitOS <- function(y, y_hat, folds, y_hat2 = NULL, level = 0.95){
  # Turning y into a matrix
  y = as.matrix(y)

  # Turning data frames into matrices
  if(is.data.frame(y_hat)){
    y_hat = as.matrix(y_hat)
  }

  if(is.data.frame(folds)){
    folds = as.matrix(folds)
  }

  if(is.data.frame(y_hat2)){
    y_hat2 = as.matrix(y_hat2)
  }

  # Sanity checks
  if(nrow(y) != nrow(y_hat) | nrow(y) != nrow(folds) | nrow(folds) != nrow(y_hat)){
    stop("y and y_hat must have the same amount of rows")
  }

  if(ncol(folds) != ncol(y_hat)){
    stop("folds and y_hat must have the same column number")
  }

  # Sanity checks for second mu
  if(!is.null(y_hat2)){
    if(ncol(folds) != ncol(y_hat2)){
      stop("folds and y_hat2 must have the same column number")
    }

    if(nrow(y) != nrow(y_hat2) | nrow(folds) != nrow(y_hat2)){
      stop("y and y_hat2 must have the same amount of rows")
    }
  }


  # number of sample splits
  n_splits <- ncol(y_hat)

  # SINGLE lapply over splits
  res_list <- lapply(seq_len(n_splits), function(i) {
    # choose appropriate y_hat2 column (or NULL)
    this_yhat2 <- if (is.null(y_hat2)) NULL else y_hat2[, i]

    # run foldOS once
    tmp <- foldOS(y, y_hat[, i], folds[, i], this_yhat2, level)

    # add split id to est
    tmp$est$split <- i

    # return everything we might need
    list(est   = tmp$est, fold1 = tmp$fold.est, fold2 = if (!is.null(y_hat2)) tmp$fold.est2 else NULL)
  })

  # helper to extract and rbind a component if it exists
  bind_component <- function(lst, name) {
    comp <- lapply(lst, `[[`, name)
    # filter out NULLs (e.g., fold2 when y_hat2 is NULL)
    comp <- comp[!vapply(comp, is.null, logical(1))]
    if (length(comp) == 0L) return(NULL)
    do.call(rbind, comp)
  }

  # Row-binding the respective list components and assigning them to variables
  est <- bind_component(res_list, "est")
  folds_one <- bind_component(res_list, "fold1")
  folds_two <- bind_component(res_list, "fold2")

  # build output
  out <- list(result = est, fold1.vals  = folds_one,n = length(y), conf.level  = level)

  # Changing output based on if yhat2 is passed
  if (!is.null(folds_two)) {
    out$fold2.vals <- folds_two
  }

  return(out)
}
