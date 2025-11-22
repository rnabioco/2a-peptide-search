#! /usr/bin/env bash

#BSUB -J hmmsearch
#BSUB -o logs/hmmsearch.out
#BSUB -e logs/hmmsearch.err
#BSUB -n 12

project="$HOME/devel/rnabioco/2a-peptide-search"
db="$project/ref/uniprot.2023_01.fa.gz"
hmm="$project/data/2a.viral.refseq.hmm"

hmmsearch --cpu 12 --max \
    -o $project/results/max/2a.uniprot.hmsearch \
    -A $project/results/max/2a.uniprot.sto \
    $hmm $db 

