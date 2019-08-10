import os
import sys

gamma_or_taus = ["tau"]
setups = ["93_95", "98_00", "urban", "rural"]
algs = ["non_const", "const", "naive"]
nreps = ["1", "2", "3","4","5","6","7","8","9","10"]


# gamma_minimax, gamma_plugin, all_tau, partial_tau_bootstrap, direct_tau_bootstrap


for alg in algs_tau:
  for setup in setups:
    for rep in nreps:
      os.system('Rscript sims_ak07.R tau ' +setup + ' ' +alg + ' 1 '+nrep)

