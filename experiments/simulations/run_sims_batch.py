import os
import sys


const_algs = ['TR', 'OLS', 'oracle_learner', 'sample_means']
const_setups = ['A','B','C','D']
non_const_algs = ['DiD_AMLE','DiD_AIPW','OLS','IPW']
non_const_setups = ['C','F','D','E']
ns = ['100','200','500','1000']
ps = ['6','12']
NREP = '100'

for alg in const_algs:
	for setup in const_setups:
		for n in ns:
			for p in ps:
				os.system('Rscript run_sims.R ' + alg + ' ' + setup + ' ' + n +
                ' ' + p + ' ' + NREP + ' constant')

for alg in non_const_algs:
    for setup in non_const_setups:
        for n in ns:
            for p in ps:
                os.system('Rscript run_sims.R ' + alg + ' ' + setup + ' ' + n +
                ' ' + p + ' ' + NREP + ' non_constant')

