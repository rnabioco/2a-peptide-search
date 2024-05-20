#! /usr/bin/env bash

#BSUB -J hmmsearch[1-2]
#BSUB -eo hmmesearch-ref-proteomes.%J.errout

set -x

class="class-$LSB_JOBINDEX"
project=$HOME/devel/rnabioco/2a-peptide-search
hmm="$project/curated-models/2A-$class.hmm"
seq=$project/data/ref/uniprot/reference_proteomes/reference_proteomes.fasta.gz

hmmsearch  \
    --tblout ref-proteomes.$class.hmmsearch.tab \
    -A ref-proteomes.$class.sto \
    --noali \
    $hmm $seq \
    | gzip -c > ref-proteomes.$class.hmmsearch.gz
