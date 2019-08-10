library(Matrix)
library(foreach)
library(magrittr)
library(grf)
library(glmnet)
library(rlearner)
library(diffindiff)
library(EQL)
library(rpart)
library(foreign)
library(sandwich)

args = commandArgs(trailingOnly = TRUE)
gamma_or_tau = as.character(args[1])
setup = as.character(args[2])
alg = as.character(args[3])
block = as.numeric(args[4])
nrep = as.numeric(args[5])

data_file = paste("./data/ak07data_",setup,".csv", sep = "")
all_data = read.csv(data_file)[,2:26]
X = all_data[,1:9]
Ti = X[,7]
Si = X[,8]
Y = X[,9]
X = X[,1:6]

# first check whether running gamma or tau
if (gamma_or_tau == "gamma"){
  if (alg == "amle"){
    # if alg is amle, this means we want to get minimax weights
    # experiment shows sherlock can handel 3000 data points at a time
    # so we split the data into blocks, and just run on a specific block
    # which is specified by the input variable block
    # 93_95 has 6, 98_00 has 10, urban has 5, rural has 2
    if (setup == "93_95"){
      num_block = 6
    } else if (setup == "98_00"){
      num_block = 10
    } else if (setup == "urban"){
      num_block = 5
    } else if (setup == "rural"){
      num_block = 2
    } else {
      stop("bad block")
    }

    if (block < num_block){
        X = X[((block-1)*3000 + 1):(block*3000),]
        Ti = Ti[((block-1)*3000 + 1):(block*3000)]
        Si = Si[((block-1)*3000 + 1):(block*3000)]
    } else if (block == num_block) {
      if (setup == "93_95"){
        X = X[15001:18614,]
        Ti = Ti[15001:18614]
        Si = Si[15001:18614]
      } else if (setup == "98_00"){
        X = X[27001:29254,]
        Ti = Ti[27001:29254]
        Si = Si[27001:29254]
      } else if (setup == "urban"){
        X = X[12001:13367,]
        Ti = Ti[12001:13367]
        Si = Si[12001:13367]
      } else if (setup == "rural"){
        X = X[3001:5247,]
        Ti = Ti[3001:5247]
        Si = Si[3001:5247]
      } else {
        stop("bad block")
      }
    }

    make_matrix = function(x) stats::model.matrix(~.-1, x)
    X = data.frame(X) %>% make_matrix
    order = 1
    dd = min(1000, nrow(X) - ncol(X))
    while (dd > 0 & order < 10) {
      order = order + 1
      dd = dd - choose(ncol(X) - 1 + order, order)
    }
    print("gotten order")
    Basis = generate.basis(X,order)
    Basis = cbind(1,Basis)
    Basis = t(Basis)
    Basis = Basis[complete.cases(Basis),]
    Basis = t(Basis)
    print("have basis")
    gamma.minimax = minimax(Basis, Ti, Si)
    filenm = paste("./real_results/results_", setup,"/gamma_minimax/gamma_", as.character(block),".csv", sep="")
    write.csv(gamma.minimax, filenm)
  } else if (alg == "aipw") {
    # in this case the gammas should be fitted
    make_matrix = function(x) stats::model.matrix(~.-1, x)
    X = data.frame(X) %>% make_matrix
    order = 1
    dd = min(1000, nrow(X) - ncol(X))
    while (dd > 0 & order < 10) {
      order = order + 1
      dd = dd - choose(ncol(X) - 1 + order, order)
    }
    print("gotten order")
    Basis = generate.basis(X,order)
    Basis = cbind(1,Basis)
    Basis = t(Basis)
    Basis = Basis[complete.cases(Basis),]
    Basis = t(Basis)
    print("have basis")
    gamma.plugin = plugin(Basis, Ti, Si)
    print("have plugin weights")
    filenm = paste("./real_results/results_",setup, "/gamma_plugin/gamma_plugin_",as.character(nrep),".csv",sep="")
    write.csv(gamma.plugin,filenm)
  } else {
    stop("bad alg")
  }
} else if (gamma_or_tau == "tau"){
  if (bootstrap == "bootstrap"){
    index = sample(1:dim(X)[1], dim(X)[1], replace = TRUE)
    X = X[index,]
    Ti = Ti[index]
    Si = Si[index]
    Y = Y[index]
    print(setup)
    write.csv(cbind(X,Ti,Si,Y), paste("./real_results/results_",
                                      setup,"/bootstrap_results/ak07data_boot_",
                                      as.character(nrep),".csv",sep=""))
  }
  # in this case we want to estimate tau
  make_matrix = function(x) stats::model.matrix(~.-1, x)
  X = data.frame(X) %>% make_matrix
  order = 1
  dd = min(1000, nrow(X) - ncol(X))
  while (dd > 0 & order < 10) {
    order = order + 1
    dd = dd - choose(ncol(X) - 1 + order, order)
  }
  print("gotten order")
  Basis = generate.basis(X,order)
  Basis = cbind(1,Basis)
  Basis = t(Basis)
  Basis = Basis[complete.cases(Basis),]
  Basis = t(Basis)
  print("have basis")

  # we modified the partial and direct learners to output the
  # non-parametric tau_hat, instead of an average value, so we can apply the gammas
  # manually
  if (alg == "non_const"){
    tau_hat = DiD(X=Basis, Y=Y, Ti=Ti, Gi=Si,
                              constant_eff="non_constant")
  } else if (alg == "const"){
    tau_hat = DiD(X=Basis, Y=Y, Ti=Ti, Gi=Si,
                              constant_eff="constant")
  } else if (alg == "naive"){
    tau_hat = naive_learner(X=X, Y=Y, Ti=Ti, Gi = Si, constant_eff = "constant")
  }


  subflnm = "/tau_results/"
  filenm = paste("./real_results/results_",setup, subflnm, alg,as.character(nrep),bootstrap,".csv", sep="")
  write.csv(tau_hat, filenm)
} else {
  stop("gamma or tau?")
}




