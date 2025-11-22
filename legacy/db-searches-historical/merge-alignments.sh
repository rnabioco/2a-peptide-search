#! /usr/bin/env bash

set -x

models=(class-1 class-2)

for model in ${models[@]}; do
    ls *$model.sto.gz > $model.txt
    esl-alimerge --list $model.txt > combined-$model.sto
    gzip combined-$model.sto
    rm -f $model.txt
done

