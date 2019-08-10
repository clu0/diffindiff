#! /bin/bash

constalgs=('TR' 'OLS' 'oracle_learner' 'sample_means')
constsetups=('A' 'B' 'C' 'D')
nonconstalgs=('DiD_AMLE' 'DiD_AIPW' 'OLS' 'IPW')
nonconstsetups=('C' 'F' 'D' 'E')
ns=(100 200 500 1000)
ps=(6 12)
rep=100

for ((i1=0; i1<${#constalgs[@]}; i1++))
do
for ((i2=0; i2<${#constsetups[@]}; i2++))
do
for ((i3=0; i3<${#ns[@]}; i3++))
do
for ((i4=0; i4<${#constalgs[@]}; i4++))
do
  alg=${constalgs[$i1]}
  setup=${constsetups[$i2]}
  n=${ns[$i3]}
  p=${ps[$i4]}

  fnm="logging/progress-$alg-$setup-$n-$p.out"
  echo $fnm

  Rscript run_sims.R $alg $setup $n $p $rep 2>&1 | tee $fnm &
done
done
done
done

for ((i1=0; i1<${#nonconstalgs[@]}; i1++))
do
for ((i2=0; i2<${#nonconstsetups[@]}; i2++))
do
for ((i3=0; i3<${#ns[@]}; i3++))
do
for ((i4=0; i4<${#constalgs[@]}; i4++))
do
  alg=${constalgs[$i1]}
  setup=${constsetups[$i2]}
  n=${ns[$i3]}
  p=${ps[$i4]}

  fnm="logging/progress-$alg-$setup-$n-$p.out"
  echo $fnm

  Rscript run_sims.R $alg $setup $n $p $rep 2>&1 | tee $fnm &
done
done
done
done
