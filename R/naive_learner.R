#' @title baseline methods used in simulations
#' @description two simple baselines: if the treatment effect is assumed to be
#' constant, then we estimate ATE by running OLS; if the treatment effect is not
#' assumed to be constant, then ATE is estimated by taking sample averages of
#' the outcomes Y (details see Lu, Nie, Wager (2019))
#'
#' @param X input features
#' @parem Y the observed responses
#' @param Ti the time variables (0 or 1)
#' @param Si the state variables (0 or 1)
#' @param constant_eff whether we are assuming a constant or non-constant
#' treatment effect (choose from "constant" and "non_constant")
#'
#'
#' @return estimated ATE and standard error
#'
#' @export

naive_learner = function(X, Y, Ti, Si, constant_eff = c("constant","non_constant")) {
  constant_eff = match.arg(constant_eff)
  if (constant_eff == "constant"){
    n = dim(X)[1]
    p = dim(X)[2]
    input_data = as.data.frame(cbind(X, Ti, Si, Ti*Si, Y))
    name_root = rep("x",p)
    names = sapply(1:p, function(i){
      name_root[i] = paste(name_root[i],toString(i),sep="")
    })
    colnames(input_data) = c(names,"Ti","Si","TiSi","Y")
    tau_fit = lm(Y~., input_data)
    tau_fit_summary = summary(tau_fit)
    TAU_hat = tau_fit_summary$coefficients[(1+p+3),1]
    std.err.est.def = tau_fit_summary$coefficients[(1+p+3),2]
    std.err.est = sqrt(sandwich::vcovHC(tau_fit)[10,10])
    tau_hat = tau_fit$coefficients[10]

    return(list(TAU_hat = tau_hat, std.err.est = std.err.est))
  } else if (constant_eff == "non_constant"){
    Y00 = mean(Y[((1-Ti)*(1-Si)==1)])
    Y00_var = var(Y[((1-Ti)*(1-Si)==1)])/sum(as.numeric((1-Ti)*(1-Si)==1))
    Y10 = mean(Y[(Ti*(1-Si)==1)])
    Y10_var = var(Y[(Ti*(1-Si)==1)])/sum(as.numeric(Ti*(1-Si)==1))
    Y01 = mean(Y[((1-Ti)*Si==1)])
    Y01_var = var(Y[((1-Ti)*Si==1)])/sum(as.numeric((1-Ti)*Si==1))
    Y11 = mean(Y[(Ti*Si==1)])
    Y11_var = var(Y[(Ti*Si==1)])/sum(as.numeric(Ti*Si==1))

    naive_tau = Y11 - Y10 - Y01 + Y00
    naive_sd = sqrt(Y00_var+Y01_var+Y10_var+Y11_var)
    return(list(TAU_hat = naive_tau, std.err.est = naive_sd))
  } else {
    stop("Effect constant of non-constant?")
  }
}
