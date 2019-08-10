#! /bin/bash

algs=('amle' 'aipw')
blocks9395=(1 2 3 4 5 6)
blocks9800=(1 2 3 4 5 6 7 8 9 10)
blocksurban=(1 2 3 4 5)
blocksrural=(1 2)

for ((i1=0; i1<${#algs[@]}; i1++))
do
  alg=${algs[$i1]}
  for ((b1=0; b1<${#blocks9395}; b1++))
  do
    block=${blocks9395[$b1]}
    fnm="logging/progress-93_95-$alg-$block-1.out"
    echo $fnm
    Rscript sims_ak07.R gamma 93_95 $alg $block 1 2>&1 | tee $fnm &
  done
  for ((b2=0; b2<${#blocks9800}; b2++))
  do
    block=${blocks9800[$b2]}
    fnm="logging/progress-98_00-$alg-$block-1.out"
    echo $fnm
    Rscript sims_ak07.R gamma 98_00 $alg $block 1 2>&1 | tee $fnm &
  done
  for ((b3=0; b3<${#blocksurban}; b3++))
  do
    block=${blocksurban[$b3]}
    fnm="logging/progress-urban-$alg-$block-1.out"
    echo $fnm
    Rscript sims_ak07.R gamma urban $alg $block 1 2>&1 | tee $fnm &
  done
  for ((b4=0; b4<${#blocksrural}; b4++))
  do
    block=${blocksrural[$b4]}
    fnm="logging/progress-rural-$alg-$block-1.out"
    echo $fnm
    Rscript sims_ak07.R gamma rural $alg $block 1 2>&1 | tee $fnm &
  done
done