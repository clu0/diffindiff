#' @title oracle methods
#' @description function that takes in oracle values for all nuisance parameters,
#' and estimates the ATE, either using the transformed regression (TR), if we assume a constant effect, or using
#' non-parametric regression (DiD-AMLE), if we assume a non-constant effect (details see Lu, Nie, Wager (2019))
#'
#' @param X input features
#' @param Y observed responses
#' @param Ti time variable
#' @param Si state variable
#' @param f the effect of Ti on outcome:
#' f(x) = E[Yi | Xi = x, Ti = 1, Si = 0] - E[Yi | Xi = x, Ti = 0, Si = 0]
#' @param h the effect of Si on outcome:
#' h(x) = E[Yi | Xi = x, Ti = 0, Si = 1] - E[Yi | Xi = x, Ti = 0, Si = 0]
#' @param b the baseline effect of input features: b(x) = E[Yi | Xi = 1, Ti = 0, Si = 0]
#' @param TAU the true ATE
#' @param tau the heterogeneous treatment effect, given as a vector with the same
#' length as Y; left as NULL if we assume the treatment is non-constant
#' @param gamma optional balancing weights
#' @param P_Ta_Sb the conditional probabilities:
#' P_Ta_Sb(x) = P(Ti = a, Si = b | Xi = x), where a and b take values in 0 and 1
#' @param constant_eff whether we assume a constant or non-constant treatment effect
#' @param lambda_choice how to choose the lasso estimate from glmnet;
#' choose from "lambda.min" and "lambda.1se"; default "lambda.min"
#'
#' @return estimated ATE
#' @export

oracle_learner = function(X, Y, Ti, Si, f, h, b, TAU, tau = NULL, gamma = NULL,
                      P_T1_S1, P_T1_S0, P_T0_S1, P_T0_S0, constant_eff = c("constant", "non_constant"), lambda_choice = "lambda.min") {
  constant_eff = match.arg(constant_eff)
  P_T1 = P_T1_S0 + P_T1_S1
  P_T0 = P_T0_S0 + P_T0_S1
  P_S1 = P_T1_S1 + P_T0_S1
  P_S0 = P_T1_S0 + P_T0_S0
  P_T1_giv_S1 = P_T1_S1/P_S1
  P_T1_giv_S0 = P_T1_S0/P_S0
  P_T0_giv_S1 = P_T0_S1/P_S1
  P_T0_giv_S0 = P_T0_S0/P_S0
  P_S1_giv_T1 = P_T1_S1/P_T1
  P_S1_giv_T0 = P_T0_S1/P_T0
  P_S0_giv_T1 = P_T1_S0/P_T1
  P_S0_giv_T0 = P_T0_S0/P_T0

  delta_T_G = P_T1_giv_S1 - P_T1_giv_S0
  delta_G_T = P_S1_giv_T1 - P_S1_giv_T0
  A  = (1/(1-delta_T_G*delta_G_T))*(Ti - P_T1 - delta_T_G*(Si - P_S1))
  B  = (1/(1-delta_T_G*delta_G_T))*(Si - P_S1 - delta_G_T*(Ti - P_T1))
  C = Ti*Si - P_T1_S1 - P_S1_giv_T1*A - P_T1_giv_S1*B

  if (is.null(tau)){
    nu = f + P_S1_giv_T1*TAU + delta_G_T*h
    sigma = h + P_T1_giv_S1*TAU + delta_T_G*f
  } else {
    nu = f + P_S1_giv_T1*tau + delta_G_T*h
    sigma = h + P_T1_giv_S1*tau + delta_T_G*f
  }

  if (is.null(tau)){
    m = b + P_T1*f + P_S1*h + P_T1_S1*TAU
  } else {
    m = b + P_T1*f + P_S1*h + P_T1_S1*tau
  }

  S = Y - m - A*nu - B*sigma

  if (constant_eff == "constant"){
    input_data = as.data.frame(cbind(S,C))
    colnames(input_data) = c("S", "C")
    tau_fit_oracle = lm(S~C - 1, input_data)

    return(list(TAU_hat = tau_fit_oracle$coefficients, std.err.est = NA))
  } else if (constant_eff == "non_constant"){
    x_tilde = cbind(as.numeric(C) * cbind(1, X))
    p.fac = rep(1,dim(x_tilde)[2])
    p.fac[1] = 0
    x_pred = cbind(1, X)
    n = length(Ti)
    k_folds = 10
    foldid = sample(rep(seq(k_folds), length = length(Ti)))
    tau_fit = glmnet::cv.glmnet(x_tilde, S, foldid = foldid, penalty.factor = p.fac)
    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    tau_hat = x_pred %*% tau_beta

    TAU_hat = mean(tau_hat + gamma*(Y - (b + Ti*f + Si*h + Ti*Si*tau_hat)))

    return(list(TAU_hat = TAU_hat, std.err.est = NA))
  } else {
    stop("what effect?")
  }
}
