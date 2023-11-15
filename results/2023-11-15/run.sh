#! /usr/bin/env bash

#BSUB -J hmmsearch[1-2]
#BSUB -eo logs/hmmsearch.%J.%I.errout
#BSUB -n 12

set -x

project=$HOME/devel/rnabioco/2a-peptide-search
hmm=$project/curated-models/2A-class-${LSB_JOBINDEX}.hmm
fa=$project/data/ref/uniprot/uniprot.2023_01.fa.gz
outbase="uniprot.class-${LSB_JOBINDEX}"

hmmsearch --tblout $outbase.tbl -o $outbase.hmmsearch -A $outbase.sto --cpu 12 $hmm $fa

gzip *.sto *.tbl *.hmmsearch

