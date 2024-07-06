#! /usr/bin/env bash

models=(2A-class-1 2A-class-2)
fa=IMGVR_all_proteins-high_confidence.fixed.faa.gz

for model in ${models[@]}; do
    hmm=../../curated-models/$model.hmm
    hmmsearch -A imgvr.$model.sto --tblout imgvr.$model.tbl.hmmsearch $hmm $fa | gzip -c > imgvr.$model.hmmsearch.gz
done
