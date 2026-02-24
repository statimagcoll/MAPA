# These functions are essentially mapaFolds and mapa but adapted to binary phenotypes.
# Can work on combining into cohesive functions that handle continuous or binary

#' @export
mapaFolds_auc <- function(y, y_hat, folds, y_hat2 = NULL, level = 0.95, ...){
  #browser()
  # Checking that variable lengths match
  if(length(y_hat) != length(y) | length(y_hat) != length(folds) | length(y) != length(folds)){
    stop("y, y_hat, and folds must be of the same length.")
  }

  # Coerce y to 0,1
  y <- as.numeric(factor(y)) - 1

  # Calculate full sample level quantities (pi_hat, density functions)
  pi_hat = mean(y)
  # conditional density estimation done on the replicate level
  dens1 <- density(y_hat[y == 1], from = 0, to = 1)
  dens0 <- density(y_hat[y == 0], from = 0, to = 1)
  f1 = approxfun(dens1$x, dens1$y, rule = 2)
  f0 = approxfun(dens0$x, dens0$y, rule = 2)
  # CDFs
  F1 = ecdf(y_hat[y == 1])
  F0 = ecdf(y_hat[y == 0])

  # Creating storage vector for results and variables
  os_folds = c()
  inf_folds_direct = c()
  inf_folds_mu = c()
  n = length(y)

  # Looping through folds
  for(fold in sort(unique(folds))){
    # Calculating and storing OS estimator
    os_res = auc_os_fold(y_hat[which(folds == fold)], y[which(folds == fold)], pi_hat = pi_hat, F0 = F0, F1 = F1, f0 = f0, f1 = f1, ...)
    os_folds = c(os_folds, os_res)

    # Calculating and storing the influence vector (both direct and indirect)
    inf_res_direct = if_auc_direct(y_hat[which(folds == fold)], y[which(folds == fold)], pi_hat = pi_hat, F0 = F0, F1 = F1, f0 = f0, f1 = f1, ...)
    inf_folds_direct = c(inf_folds_direct, inf_res_direct)
    inf_res_mu = if_auc_mu(y_hat[which(folds == fold)], y[which(folds == fold)], f0 = f0, f1 = f1)
    inf_folds_mu = c(inf_folds_mu, inf_res_mu)
  }

  # Calculating the mean OS estimate across folds
  os_estimate = mean(os_folds, na.rm = T)

  # Calculating the standard error of the OS estimate
  stand_error = sqrt(var(inf_folds_direct + inf_folds_mu)/n)

  # Calculating the bounds of the OS estimate
  LB_est = os_estimate - qnorm((1 + level) / 2)*stand_error
  UB_est = os_estimate - qnorm((1 - level) / 2)*stand_error

  # Denoting the fold number that each OS fold estimate belongs to
  names(os_folds) <- sort(unique(folds))

  # Creating a results data frame for output
  estimate = data.frame(est = os_estimate, LB = LB_est, UB = UB_est, se = stand_error)


  if(!is.null(y_hat2)){
    # Checking that variable lengths match
    if(length(y_hat2) != length(y) | length(y_hat2) != length(folds) | length(y) != length(folds)){
      stop("y, y_hat2, and folds must be of the same length.")
    }

    # conditional density estimation done on the replicate level
    dens1_2 <- density(y_hat2[y == 1], from = 0, to = 1)
    dens0_2 <- density(y_hat2[y == 0], from = 0, to = 1)
    f1_2 = approxfun(dens1$x, dens1$y, rule = 2)
    f0_2 = approxfun(dens0$x, dens0$y, rule = 2)
    # CDFs
    F1_2 = ecdf(y_hat2[y == 1])
    F0_2 = ecdf(y_hat2[y == 0])

    # Creating storage vector for results and variables
    os_folds2 = c()
    inf_folds2_direct = c()
    inf_folds2_mu = c()
    n = length(y)

    # Looping through folds
    for(fold in sort(unique(folds))){
      # Calculating and storing OS estimator
      os_res = auc_os_fold(y_hat2[which(folds == fold)], y[which(folds == fold)], pi_hat = pi_hat, F0 = F0_2, F1 = F1_2, f0 = f0_2, f1 = f1_2, ...)
      os_folds2 = c(os_folds2, os_res)

      # Calculating and storing the influence vector (both direct and indirect)
      inf_res_direct = if_auc_direct(y_hat2[which(folds == fold)], y[which(folds == fold)], pi_hat = pi_hat, F0 = F0_2, F1 = F1_2, f0 = f0_2, f1 = f1_2, ...)
      inf_folds2_direct = c(inf_folds2_direct, inf_res_direct)
      inf_res_mu = if_auc_mu(y_hat2[which(folds == fold)], y[which(folds == fold)], f0 = f0_2, f1 = f1_2)
      inf_folds2_mu = c(inf_folds2_mu, inf_res_mu)
    }

    # Calculating the mean OS estimate across folds
    os_estimate2 = mean(os_folds2, na.rm = T)

    # Calculating standard error
    stand_error2 =  sqrt(var(inf_folds2_direct + inf_folds2_mu)/n)

    # Calculating the OS estimate bounds across folds
    LB_est2 = os_estimate2 - qnorm((1 + level) / 2)* stand_error2
    UB_est2 = os_estimate2 - qnorm((1 - level) / 2)*stand_error2

    # Denoting the fold number that each OS fold estimate belongs to
    names(os_folds2) <- sort(unique(folds))

    # Creating data frame output
    estimate2 = data.frame(est = os_estimate2, LB = LB_est2, UB = UB_est2, se = stand_error2)

    # Calculating OS difference
    diff = auc_os_diff(os_estimate2, os_estimate, inf_folds_direct + inf_folds_mu, inf_folds2_direct + inf_folds2_mu)
    os_diff = data.frame(est = diff[1], LB = diff[2], UB = diff[3], se = diff[4])

    # Combining all data into one
    est_total = rbind(estimate, estimate2, os_diff)

    # Assigning row names
    est_total$metric <- c("yhat1", "yhat2", "Diff")
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

#' @export
mapa_auc <- function(y, y_hat, folds, y_hat2 = NULL, level = 0.95, ...){
  #browser()
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


  # Assigning the number of sample splits
  n_splits <- ncol(y_hat)


  # Gathering results
  res_list <- lapply(seq_len(n_splits), function(i) {
    tmp <- mapaFolds_auc(y, y_hat[, i], folds[, i], y_hat2[,i], level, ...)
    tmp$est$split <- i

    list(
      est   = tmp$est,
      fold1 = tmp$fold.est,
      fold2 = if (is.null(y_hat2)) NULL else tmp$fold.est2
    )
  })

  # Combining all results
  out <- list(
    est        = do.call(rbind, lapply(res_list, `[[`, "est")),
    est1.folds = do.call(rbind, lapply(res_list, `[[`, "fold1")),
    n          = length(y),
    conf.level = level
  )

  # Naming rows in est1.folds output to signal sample split number
  rownames(out$est1.folds) <- paste0("ss", seq_len(nrow(out$est1.folds)))

  if (!is.null(y_hat2)) {
    # Adding est2.folds when yhat2 is present
    out$est2.folds <- do.call(rbind, lapply(res_list, `[[`, "fold2"))

    # Naming rows in est2.folds output to signal sample split number
    rownames(out$est2.folds) <- paste0("ss", seq_len(nrow(out$est2.folds)))
  }

  # Adding a summary data frame element in case there is more than one sample split
  out$est_summary = out$est %>%
    group_by(metric) %>%
    summarise(
      est = median(est, na.rm = TRUE),
      LB = median(LB, na.rm = TRUE),
      UB = median(UB, na.rm = TRUE),
      se = median(se, na.rm = TRUE)
    )

  # Storing est_summary as a data frame instead of a tibble
  out$est_summary = data.frame(out$est_summary)

  # Returning output
  return(out)
}
