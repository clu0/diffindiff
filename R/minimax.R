#' @title Minimax balancing weights
#' @description function that finds the minimax balancing weights, using the method in
#' Lu, Nie, Wager (2019)
#'
#' @param X the input features
#' @param Ti the time variable (0 or 1)
#' @param Si the state variable (0 or 1)
#' @param zeta tuning parameter, in [0,1]
#' @param solver cvx opmitization solver, choose from "ECOS" or "SCS"
#'
#' @return a vector gamma of the same length as Ti
#'
#' @export


minimax = function(X, Ti, Si, zeta=0.5, solver = c("ECOS", "SCS"), verbose = FALSE) {
  solver = match.arg(solver)
  nobs = nrow(X)
  pobs = ncol(X)
  gg = CVXR::Variable(nobs + 4)
  objective = (1 - zeta) * sum(gg[1:nobs]^2) + zeta * sum(gg[nobs + 1:4]^2)
  contraints = list(
    sum(gg[1:nobs]) == 0,
    t(X) %*% gg[1:nobs] <= gg[nobs + 1],
    -t(X) %*% gg[1:nobs] <= gg[nobs + 1],
    t(X) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
    -t(X) %*% (Ti * gg[1:nobs]) <= gg[nobs + 2],
    t(X) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
    -t(X) %*% (Si * gg[1:nobs]) <= gg[nobs + 3],
    sum(Ti * Si * gg[1:nobs]) == 1,
    t(X) %*% (Ti * Si * gg[1:nobs]) <= colMeans(X) + gg[nobs + 4],
    - t(X) %*% (Ti * Si * gg[1:nobs]) <= - colMeans(X) + gg[nobs + 4]
  )
  cvx.problem = CVXR::Problem(CVXR::Minimize(objective), contraints)
  cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
  result = cvx.output$getValue(gg)
  gamma = nobs * result[1:nobs]
  return(gamma)
}
