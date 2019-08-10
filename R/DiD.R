#' @title non-param diff-in-diff
#'
#' @description  Non-parametric diff-in-diff estimation, as proposed by Lu, Nie, Wager (2019),
#'  implemented with the lasso (glmnet) and the R learner (rlearner) from Nie and Wager (2017)
#'
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param Ti the time variable (0 or 1)
#' @param Si the state variable (0 or 1)
#' @param constant_eff whether the treatment effect is modeled to be constant, choose from "constant" or "non_constant"
#' if assuming a constant treatment effect, set constant_eff = "constant" (which runs the Transformed Regression from the paper);
#' if assuming a non-constant effect, set constant_eff = "non_constant"
#' @param gamma balancing weights, a vector of the same length as the length of Y (observables),
#' used for balancing in the treatment estimation; can be left blank, in which case the algorithm
#' does not use balancing
#' @param k_folds number of folds for cross-validation; default is 10
#' @param lambda_choice choose how to cross-validate; choose from "lambda.min" or "lambda.1se"
#' @param penalty_factor user specified penalty factor, used in the R learner
#'
#' @return a list including TAU_hat, estimated ATE; std.err.est, estimated standard error;
#' tau_hat, non-parametric point estimates of treatment effect; y.est, estimated outcome values.
#'
#' @export

DiD = function(X, Y, Ti, Si, constant_eff = c("constant","non_constant"), gamma = NULL, k_folds = 10,
                           lambda_choice = "lambda.min", penalty_factor = NULL) {
  constant_eff = match.arg(constant_eff)
  n = length(Ti)
  if (is.null(penalty_factor)){
    penalty_factor = rep(1, dim(X)[2])
  }

  # create foldid
  foldid = sample(rep(seq(k_folds), length = length(Ti)))

  # Fit values
  nu_fit = rlearner::rlasso(X, Ti, Y, k_folds = k_folds, penalty_factor = penalty_factor)
  nu_hat = nu_fit$tau_hat
  t_hat = nu_fit$p_hat
  m_hat = nu_fit$m_hat

  sigma_fit = rlearner::rlasso(X, Si, Y, k_folds = k_folds, penalty_factor = penalty_factor)
  sigma_hat = sigma_fit$tau_hat
  s_hat = sigma_fit$p_hat

  delta_fit = glmnet::cv.glmnet(X, (Ti - t_hat) * (Si - s_hat), foldid = foldid, keep = TRUE, alpha = 1)
  delta_hat = as.matrix(delta_fit$fit.preval[,!is.na(colSums(delta_fit$fit.preval))])[, delta_fit$lambda == delta_fit$lambda.min]

  scaling = 1 - (delta_hat^2 / (s_hat *(1 - s_hat) * t_hat * (1 - t_hat)))
  A_hat = (Ti - t_hat - (delta_hat * (Si - s_hat)) / (s_hat * (Si - s_hat))) / scaling
  B_hat = (Si - s_hat - (delta_hat * (Ti - t_hat)) / (t_hat * (Ti - t_hat))) / scaling
  C_hat = Si * Ti - (s_hat + delta_hat / t_hat) * A_hat - (t_hat + delta_hat / s_hat) * B_hat - (s_hat * t_hat + delta_hat)

  S_hat = Y - m_hat - A_hat*nu_hat - B_hat*sigma_hat

  if (constant_eff == "constant") {
    input_data = as.data.frame(cbind(S_hat,C_hat))
    colnames(input_data) = c("S_hat", "C_hat")
    tau_fit_model = lm(S_hat~C_hat - 1, input_data)
    tau_fit_summary = summary(tau_fit_model)
    TAU_hat = tau_fit_summary$coefficients[1]
    std.err.est = sqrt(sandwich::vcovHC(tau_fit_model))

    return(list(TAU_hat = TAU_hat, std.err.est = std.err.est))
  } else if (constant_eff == "non_constant"){
    x_tilde = cbind(as.numeric(C_hat) * cbind(1, X))
    p.fac = rep(1,dim(x_tilde)[2])
    p.fac[1] = 0
    x_pred = cbind(1, X)
    n = length(Ti)
    k_folds = 10
    foldid = sample(rep(seq(k_folds), length = length(Ti)))
    tau_fit = glmnet::cv.glmnet(x_tilde, S_hat, foldid = foldid, penalty.factor = p.fac)
    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    tau_hat = x_pred %*% tau_beta

    if (is.null(gamma)){
      TAU_hat = mean(tau_hat)
      std.err.est = NA
    } else {
      y.est = m_hat + A_hat*nu_hat + B_hat*sigma_hat + C_hat*tau_hat
      TAU_hat = mean(tau_hat + gamma*(Y - y.est))
      std.err.est = mean((tau_hat - TAU_hat)^2 + gamma^2 * (Y - y.est)^2)
    }
    y.est = m_hat + A_hat*nu_hat + B_hat*sigma_hat + C_hat*tau_hat
    return(list(TAU_hat = TAU_hat, std.err.est = std.err.est,
                tau_hat = tau_hat, y.est = y.est))
  } else {
    stop("what effect?")
  }
}
