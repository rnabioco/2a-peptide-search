#! /usr/bin/env bash

models=(2A-class-1 2A-class-2)
fa=new-old.fa

for model in ${models[@]}; do
    hmm=../../curated-models/$model.hmm
    hmmsearch -A new-old.$model.sto --tblout new-old.$model.tbl.hmmsearch $hmm $fa | gzip -c > new-old.$model.hmmsearch.gz
done
