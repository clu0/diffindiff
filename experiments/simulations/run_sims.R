# Simulations for diff and diff
library(Matrix)
library(foreach)
library(magrittr)
library(grf)
library(glmnet)
library(rlearner)
library(diffindiff)
#library(amlinear)
library(EQL)
library(rpart)
source("utils.R")

start_time <- Sys.time()

args = commandArgs(trailingOnly = TRUE)
alg = as.character(args[1])
setup = as.character(args[2])
n = as.numeric(args[3])
p = as.numeric(args[4])
NREP = as.numeric(args[5])
constant_eff = as.character(args[6])


eta = 0.1

if (setup == 'A') {
  get.params = function(){
    X = matrix(rnorm(n * p), n, p)
    TAU = 1
    f = X[,3]
    h = X[,4]
    S.treat.prob = 0.4
    T.treat.prob = 0.4
    Si = rbinom(n, 1, S.treat.prob)
    Ti = rbinom(n, 1, T.treat.prob)
    P_T1_S1 = T.treat.prob*S.treat.prob
    P_T1_S0 = T.treat.prob*(1-S.treat.prob)
    P_T0_S1 = (1-T.treat.prob)*S.treat.prob
    P_T0_S0 = (1-T.treat.prob)*(1-S.treat.prob)
    b = pmax(X[,1] + X[,2], 0)
    Y = b + Ti*f + Si*h + Ti*Si*TAU + rnorm(n)
    list(X=X, TAU=TAU, tau=TAU, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }

} else if (setup == 'B'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    TAU = 1
    f = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    h = log(1 + exp(-X[,4] - X[,5]))
    probabilities = matrix(0,4,n)
    locations = matrix(0,4,n)
    for (i in 1:n) {
      T1_S1 = pmax(eta, 0.7*(1 / (1+exp(-X[i,3])) ) )
      T1_S0 = pmax(eta, (0.8-T1_S1)*(1 / (1 + exp(-X[i,4]))) )
      T0_S1 = pmax(eta, (0.9-T1_S1-T1_S0)*(1 / (1 + exp(-X[i,5]))) )
      T0_S0 = 1 - T1_S1 - T1_S0 - T0_S1
      probabilities[,i] = c(T1_S1, T1_S0, T0_S1, T0_S0)
      locations[,i] = rmultinom(1, 1, c(T1_S1, T1_S0, T0_S1, T0_S0))
    }
    P_T1_S1 = probabilities[1,]
    P_T1_S0 = probabilities[2,]
    P_T0_S1 = probabilities[3,]
    P_T0_S0 = probabilities[4,]
    Ti = locations[1,] + locations[2,]
    Si = locations[1,] + locations[3,]
    b = pmax(X[,1] + X[,2], 0)
    Y = b + Ti*f + Si*h + Ti*Si*TAU + rnorm(n)
    list(X=X, TAU=TAU, tau=TAU, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else if (setup == 'E'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    tau = (X[,1] + X[,2])^2
    TAU = 2
    f = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    h = log(1 + exp(-X[,4] - X[,5]))
    probabilities = matrix(0,4,n)
    locations = matrix(0,4,n)
    for (i in 1:n) {
      T1_S1 = pmax(eta, 0.7*(1 / (1+exp(-X[i,3])) ) )
      T1_S0 = pmax(eta, (0.8-T1_S1)*(1 / (1 + exp(-X[i,4]))) )
      T0_S1 = pmax(eta, (0.9-T1_S1-T1_S0)*(1 / (1 + exp(-X[i,5]))) )
      T0_S0 = 1 - T1_S1 - T1_S0 - T0_S1
      probabilities[,i] = c(T1_S1, T1_S0, T0_S1, T0_S0)
      locations[,i] = rmultinom(1, 1, c(T1_S1, T1_S0, T0_S1, T0_S0))
    }
    P_T1_S1 = probabilities[1,]
    P_T1_S0 = probabilities[2,]
    P_T0_S1 = probabilities[3,]
    P_T0_S0 = probabilities[4,]
    Ti = locations[1,] + locations[2,]
    Si = locations[1,] + locations[3,]
    b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3])) + sin(pi*X[,1]*X[,2])
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else if (setup == 'F'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    tau = sin(2*pi*X[,1]) + X[,4] + 0.5*X[,5]
    TAU = 0
    f = 1/(1+exp(X[,3]))
    h = 1/(1+exp(X[,4]))
    probabilities = matrix(0,4,n)
    P_T1 = 0.45
    P_S1 = pmin(pmax(eta, (1 / (1+exp(-0.5*X[,3] + X[,5])))), 1-eta)
    P_T1_S1 = P_T1*P_S1
    P_T1_S0 = P_T1*(1-P_S1)
    P_T0_S1 = (1-P_T1)*P_S1
    P_T0_S0 = (1-P_T1)*(1-P_S1)
    Ti = rbinom(n,1,P_T1)
    Si = rbinom(n,1,P_S1)
    b = pmax(X[,1] + X[,2], 0)
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else if (setup == 'D'){
  get.params = function(){
    X = matrix(rnorm(n * p), n, p)
    tau = 2
    TAU = 2
    f = 5 * (sin(pi*X[,1]*X[,2]) + 2*(X[,3]-0.5)^2 + X[,4] + 0.5*X[,5])
    h = 0
    probabilities = matrix(0,4,n)
    P_T1 = 0.5
    P_S1 = 0.5
    P_T1_S1 = P_T1*P_S1
    P_T1_S0 = P_T1*(1-P_S1)
    P_T0_S1 = (1-P_T1)*P_S1
    P_T0_S0 = (1-P_T1)*(1-P_S1)
    Ti = rbinom(n,1,P_T1)
    Si = rbinom(n,1,P_S1)
    b = 0
    Y = b + Ti*f + Si*h + Ti*Si*tau + rnorm(n)
    return(list(X=X, TAU=TAU, tau=tau, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0))
  }
} else if (setup == 'C'){
  get.params = function() {
    X = matrix(rnorm(n * p), n, p)
    TAU = 1
    f = 0
    h = 0
    probabilities = matrix(0,4,n)
    locations = matrix(0,4,n)
    for (i in 1:n) {
      T1_S1 = 0.5 + 0.5*(1-6*eta)*sin(X[i,1]*1.5)
      T1_S0 = (1 - T1_S1)/3
      T0_S1 = (1-T1_S1)/3
      T0_S0 = 1 - T1_S1 - T1_S0 - T0_S1
      probabilities[,i] = c(T1_S1, T1_S0, T0_S1, T0_S0)
      locations[,i] = rmultinom(1, 1, c(T1_S1, T1_S0, T0_S1, T0_S0))
    }
    P_T1_S1 = probabilities[1,]
    P_T1_S0 = probabilities[2,]
    P_T0_S1 = probabilities[3,]
    P_T0_S0 = probabilities[4,]
    Ti = locations[1,] + locations[2,]
    Si = locations[1,] + locations[3,]
    b = 2*sin(X[,1]*1.5)
    Y = b + Ti*f + Si*h + Ti*Si*TAU + rnorm(n)
    list(X=X, TAU=TAU, tau=TAU, f=f, h=h, b=b, Ti=Ti, Si=Si, Y=Y, P_T1_S1=P_T1_S1, P_T1_S0=P_T1_S0, P_T0_S1=P_T0_S1, P_T0_S0=P_T0_S0)
  }
} else {
  stop("bad setup")
}


err_results = matrix(0,1,NREP)
std_err_results = matrix(0,1,NREP)


for (iter in 1:NREP) {
  params.train = get.params()
  make_matrix = function(x) stats::model.matrix(~.-1, x)
  X = data.frame(params.train$X) %>% make_matrix
  order = 1
  dd = min(1000, nrow(X) - ncol(X))
  while (dd > 0 & order < 10) {
    order = order + 1
    dd = dd - choose(ncol(X) - 1 + order, order)
  }
  print("gotten order")
  Basis = generate.basis(X,order)
  #Basis = cbind(1,Basis)
  Basis = t(Basis)
  Basis = Basis[complete.cases(Basis),]
  Basis = t(Basis)
  print("have basis")
  if (constant_eff == "constant") {
    if (alg == "oracle_learner"){
      fit = oracle_learner(P_T1_S1=params.train$P_T1_S1,
                       P_T1_S0=params.train$P_T1_S0,
                       P_T0_S1=params.train$P_T0_S1,
                       P_T0_S0=params.train$P_T0_S0,
                       Ti=params.train$Ti,
                       Si=params.train$Si,
                       f=params.train$f,
                       h=params.train$h,
                       b=params.train$b,
                       TAU=params.train$TAU,
                       Y=params.train$Y,
                       constant_eff = constant_eff)
    }
    else if (alg == "OLS") {
      fit = naive_learner(X = X,
                           Y = params.train$Y,
                           Ti = params.train$Ti,
                           Si = params.train$Si,
                           constant_eff = "constant")
    } else if (alg == "TR") {
      fit = DiD(X = Basis,
                             Y = params.train$Y,
                             Ti = params.train$Ti,
                             Si = params.train$Si,
                             constant_eff = constant_eff)
    } else if (alg == "sample_means"){
        fit = naive_learner(X = X,
                            Y = params.train$Y,
                            Ti = params.train$Ti,
                            Si = params.train$Si,
                            constant_eff = "non_constant")
    }
    else {
      stop("alg not yet supported")
    }
  } else if (constant_eff == "non_constant") {
    X = params.train$X
    Y = params.train$Y
    Ti = params.train$Ti
    Si = params.train$Si
    P_T1_S1 = params.train$P_T1_S1
    P_T1_S0 = params.train$P_T1_S0
    P_T0_S1 = params.train$P_T0_S1
    P_T0_S0 = params.train$P_T0_S0
    f=params.train$f
    h=params.train$h
    b=params.train$b
    TAU = params.train$TAU
    tau = params.train$tau
    if (alg == "OLS"){
      fit = naive_learner(X = X,
                           Y = params.train$Y,
                           Ti = params.train$Ti,
                           Si = params.train$Si,
                           constant_eff = "constant")
    } else{
    fit = diffindiff.comparison(X, Y, Ti, Si, f, h, b, TAU, tau,
                         P_T1_S1, P_T1_S0, P_T0_S1, P_T0_S0, alg = alg)
    }
  }
  TAU_hat = fit$TAU_hat
  std.err.est = fit$std.err.est
  print(TAU_hat)
  print(std.err.est)
  err_results[iter] = TAU_hat - params.train$TAU
  std_err_results[iter] = std.err.est
}

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

if (constant_eff == "constant"){
  fnm = paste("results/results_const/output", alg, setup, n, p, constant_eff, NREP, "full.csv", sep="-")
} else {
  fnm = paste("results/results_non_const/output", alg, setup, n, p, constant_eff, NREP, "full.csv", sep="-")
}
results = rbind(err_results, std_err_results)
rownames(results) = c("est_err", "est_std_err")
write.csv(results, file=fnm)
print("successsfully saved")
