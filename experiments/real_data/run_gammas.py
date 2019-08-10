import os
import sys

algs = ["amle","aipw"]
blocks_93_95 = ["1", "2", "3","4","5","6"]
blocks_98_00 = ["1", "2", "3","4","5","6","7","8","9","10"]
blocks_urban = ["1", "2", "3","4","5"]
blocks_rural = ["1", "2"]

for alg in algs:
  for block in blocks_93_95:
    os.system('Rscript sims_ak07.R gamma 93_95 ' + alg + ' ' + block + ' 1')
  for block in blocks_98_00:
    os.system('Rscript sims_ak07.R gamma 98_00 ' + alg + ' ' + block + ' 1')
  for block in blocks_urban:
    os.system('Rscript sims_ak07.R gamma urban ' + alg + ' ' + block + ' 1')
  for block in blocks_rural:
    os.system('Rscript sims_ak07.R gamma rural ' + alg + ' ' + block + ' 1')

