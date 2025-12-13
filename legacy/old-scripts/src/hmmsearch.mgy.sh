#! /usr/bin/env bash

#BSUB -J hmmsearch[1-30]
#BSUB -o logs/hmmsearch.out
#BSUB -e logs/hmmsearch.err
#BSUB -n 6

project="$HOME/devel/rnabioco/2a-peptide-search"
results=$project/2023-04-14
db="$project/ref/mgy"
hmm="$project/data/2a.uniprot.hmm"

hmmsearch --cpu 6 \
    -o $project/results/2a.uniprot.hmsearch \
    -A $project/results/2a.uniprot.sto \
    $hmm $db 

