library(Matrix)
library(foreach)
library(magrittr)
library(glmnet)
library(diffindiff)
library(EQL)
library(rpart)
library(foreign)
library(latex2exp)

#args = commandArgs(trailingOnly = TRUE)
#setup = as.character(args[1])
setup = "rural"

gamma_results_header_fl = "./results/real_results/results_"
results_header_fl = "./results/real_results/results_"
data_header_fl = "./backup_real_results/ak07data_"
gamma_minimax_flnm = "/gamma_minimax/gamma_"
gamma_plugin_flnm = "/gamma_plugin"
tau_results_flnm = "/tau_results"
bootstrap_flnm = "/bootstrap_results_direct"

nrep = 5
bootstrap_nrep = 100
g_blk = 3000

get_tau_sd = function(flnm, gamma,n) {
  results = read.csv(flnm)
  TAU_hat = sum(drop(results[2]) +
                          gamma*(drop(results[3])
                                         - drop(results[4])))/n
  var.est = sum((drop(results[2])-TAU_hat)^2 +
                         (gamma^2)*
                         (drop(results[3])- drop(results[4]))^2)/n
  std.err.est = sqrt(var.est/n)
  return(list(TAU_hat = TAU_hat, std.err.est = std.err.est))
}

# get_cluster_taus = function(X, Ti, Si, Y, group_info, which_alg,
#                             tau_hat = NULL, y.est = NULL, gamma = NULL)

all_data = read.csv(paste(data_header_fl,setup,".csv",sep=""))[,2:10]#26]
X = all_data[,1:9]
Ti = X[,7]
Si = X[,8]
Y = X[,9]
X = X[,1:6]

naive_means = naive_learner(X, Y, Ti, Si, "non_constant")


n = dim(X)[1]
num_gamma = (n - n%%g_blk)/g_blk
if (setup == "93_95"){
  num_gamma = num_gamma - 1
}

gamma.minimax = rep(0,n)
for (i in 1:num_gamma){
  filenm = paste(gamma_results_header_fl,
                 setup,gamma_minimax_flnm,as.character(i),".csv", sep="")
  gamma = read.csv(filenm)[2]
  for (j in 1:g_blk){
    gamma.minimax[(i-1)*g_blk+j] = gamma[j,]
  }
}
gamma = read.csv(paste(gamma_results_header_fl,
                       setup,gamma_minimax_flnm,
                       as.numeric(num_gamma+1),".csv",sep=""))
for (j in 1:(n-g_blk*num_gamma)){
  gamma.minimax[g_blk*num_gamma+j] = gamma[j,2]
}

gamma.plugin = rep(0,n)
for (i in 1:num_gamma){
  filenm = paste(results_header_fl,
                 setup,gamma_plugin_flnm,as.character(i),".csv", sep="")
  gamma = read.csv(filenm)[2]
  for (j in 1:g_blk){
    gamma.plugin[(i-1)*g_blk+j] = gamma[j,]
  }
}
gamma = read.csv(paste(results_header_fl,
                       setup,gamma_plugin_flnm,
                       as.numeric(num_gamma+1),".csv",sep=""))
for (j in 1:(n-g_blk*num_gamma)){
  gamma.plugin[g_blk*num_gamma+j] = gamma[j,2]
}

# gamma.plugin = rep(0,n)
# for (i in 1:nrep){
#   gamma.plugin = gamma.plugin + read.csv(paste(results_header_fl,
#                                                setup,gamma_plugin_flnm,
#                                                as.character(i),".csv",sep=""))[,2]
# }
# gamma.plugin = gamma.plugin/nrep
# gamma.plugin2 = read.csv(paste(results_header_fl,setup,
#                                gamma_plugin_flnm,"1",
#                               ".csv",sep=""))[,2]

diffindiff_amles = matrix(0,2,nrep)
diffindiff_aipws = matrix(0,2,nrep)
TRs = matrix(0,2,nrep)

for (i in 1:nrep){
  diffindiff_amle = get_tau_sd(paste(results_header_fl,
                                setup, tau_results_flnm,
                                "/non_const",
                                as.character(i), "std.csv",sep=""), gamma.minimax,n)
  diffindiff_amles[1,i] = diffindiff_amle$TAU_hat
  diffindiff_amles[2,i] = diffindiff_amle$std.err.est

  diffindiff_aipw = get_tau_sd(paste(results_header_fl,
                                  setup,tau_results_flnm, "/non_const",
                                  as.character(i), ".csv",sep=""), gamma.plugin,n)
  diffindiff_aipws[1,i] = diffindiff_aipw$TAU_hat
  diffindiff_aipws[2,i] = diffindiff_aipw$std.err.est

  TR = read.csv(paste(results_header_fl,
                             setup,tau_results_flnm, "/const",
                             as.character(i), ".csv",sep=""))
  TRs[1,i] = TR[2]$TAU_hat
  TRs[2,i] = TR[3]$std.err.est
  #
  # direct_amle = get_tau_sd(paste(results_header_fl,
  #                                 setup,tau_results_flnm, "/direct_non_const",
  #                                 as.character(i), ".csv",sep=""), gamma.minimax,n)
  # direct_amles[1,i] = direct_amle$TAU_hat
  # direct_amles[2,i] = direct_amle$std.err.est
  #
  # direct_aipw = get_tau_sd(paste(results_header_fl,
  #                                 setup,tau_results_flnm, "/direct_non_const",
  #                                 as.character(i), ".csv",sep=""), gamma.plugin,n)
  # direct_aipws[1,i] = direct_aipw$TAU_hat
  # direct_aipws[2,i] = direct_aipw$std.err.est
  #
  # direct_c = read.csv(paste(results_header_fl,
  #                            setup,tau_results_flnm, "/direct_const",
  #                            as.character(i), ".csv",sep=""))
  # direct_cs[1,i] = direct_c[2]$TAU_hat
  # direct_cs[2,i] = direct_c[3]$std.err.est
}



# exploratory
partial_results = read.csv("./experiments/old_real_results/results_93_95/tau_results/partial_non_const2.csv")

hist(partial_results[,2],xlab=TeX('$\\hat{\\tau}(\\cdot)$'),xlim=c(0,0.8),main="")
# #plot(X[,4],partial_results[,2])
#
# partial_against_age = matrix(0,39,2)
# for (i in 21:59){
#   relevant_taus = partial_results[X[,6]==(i),2]
#   partial_against_age[i-20,1] = mean(relevant_taus)
#   partial_against_age[i-20,2] = sd(relevant_taus)
# }
#
# plot(21:59, partial_against_age[,1],
#      xlab="Age", ylab="Average Treatmen Effect",type="l")

partial_against_age = matrix(0, bootstrap_nrep, 39)
for (rep in 1:bootstrap_nrep){
  partial_tau_hat = read.csv(paste(results_header_fl, setup,
                                   bootstrap_flnm, "/direct_non_const",
                                   as.character(rep),"bootstrap.csv",sep=""))
  all_data = read.csv(paste(results_header_fl, setup,
                            "/bootstrap_results/ak07data_boot_",
                            as.character(rep),".csv",sep=""))[,2:10]
  X = all_data[,1:9]
  Ti = X[,7]
  Si = X[,8]
  Y = X[,9]
  X = X[,1:6]
  for (i in 21:59){
    relevant_taus = partial_tau_hat[X[,6]==(i),2]
    partial_against_age[rep,i-20] = mean(relevant_taus)
  }
}

plot(21:59, partial_against_age[1,],
     xlab="Age", ylab="Average Treatmen Effect",type="l")#,ylim=c(-0.5,3))

for (i in 2:100){
  lines(21:59, partial_against_age[i,])
}

# rural_plot = list(tau=partial_results[,2], is_rural = X[,5])
# boxplot(tau~is_rural, rural_plot, names=c("Urban", "Rural"), ylab="Tau hat")
#
# gender_plot = list(tau=partial_results[,2], is_rural = X[,2])
# boxplot(tau~is_rural, gender_plot, names=c("Female", "Male"), ylab="Tau hat")
#
# marriage_plot = list(tau=partial_results[,2], is_rural = X[,3])
# boxplot(tau~is_rural, marriage_plot, names=c("Married", "Single"))
#
# immig_plot = list(tau=partial_results[,2], is_rural = X[,4])
# boxplot(tau~is_rural, immig_plot, names=c("Non-immigrant", "Immigrant"))
#

