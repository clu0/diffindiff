#' @title Plugin weights
#' @description finds balancing weights using the plugin estimator for P(Ti=1 | X) etc;
#' plugin weights estimated using lasso (glmnet); weights are used for methods such as
#' IPW and DiD-AIPW
#'
#' @param X input features
#' @param Ti time variables
#' @param Si state variables
#' @param k_folds number of folds for cross-fitting
#'
#' @return a vector gamma.plugin of weights of the same length as Ti
#'
#' @export

plugin = function(X, Ti, Si, k_folds = 10){
  # create plug in weights:
  n = length(Ti)
  foldid = sample(rep(seq(k_folds), length = n))
  t_fit = glmnet::cv.glmnet(X, Ti, foldid = foldid, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = 1)
  t_hat = as.matrix(t_fit$fit.preval[,!is.na(colSums(t_fit$fit.preval))])[, t_fit$lambda == t_fit$lambda.min]

  s_fit = glmnet::cv.glmnet(X, Si, foldid = foldid, keep = TRUE, family = "binomial", type.measure = "deviance", alpha = 1)
  s_hat = as.matrix(s_fit$fit.preval[,!is.na(colSums(s_fit$fit.preval))])[, s_fit$lambda == s_fit$lambda.min]

  delta_fit = glmnet::cv.glmnet(X, (Ti - t_hat) * (Si - s_hat), foldid = foldid, keep = TRUE, alpha = 1)
  delta_hat = as.matrix(delta_fit$fit.preval[,!is.na(colSums(delta_fit$fit.preval))])[, delta_fit$lambda == delta_fit$lambda.min]

  P_T1_S1_hat = delta_hat + t_hat * s_hat
  P_T0_S1_hat = s_hat - P_T1_S1_hat
  P_T1_S0_hat = t_hat - P_T1_S1_hat
  P_T0_S0_hat = 1- P_T1_S1_hat - P_T1_S0_hat - P_T0_S1_hat

  gamma.plugin = (Ti*Si)/P_T1_S1_hat + (Ti*(1-Si))/P_T1_S0_hat + ((1-Ti)*Si)/P_T0_S1_hat + ((1-Ti)*(1-Si))/P_T0_S0_hat

  return(gamma.plugin)
}
