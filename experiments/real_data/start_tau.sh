#! /bin/bash

setups=('93_95' '98_00' 'urban' 'rural')
algs=('non_const' 'const' 'naive')
nreps=(1 2 3 4 5 6 7 8 9 10)

for ((i1=0; i1<${#setups[@]}; i1++))
do
for ((i2=0; i2<${algs[@]}; i2++))
do
for ((i3=0; i3<${nreps[@]}; i3++))
do

  setup=${setups[$i1]}
  alg=${algs[$i2]}
  nrep=${nreps[$i3]}

  fnm="logging/progress-$setup-$alg-1-$nrep.out"
  echo $fnm
  Rscript sims_ak07.R tau $setup $alg 1 $nrep 2>&1 | tee $fnm &
done
done
done