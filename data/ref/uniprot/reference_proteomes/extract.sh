#! /usr/bin/env bash

#BSUB -J extract-ref
#BSUB -eo logs/extract-ref.%J.errout

# tar zxvf Reference_Proteomes_2024_02.tar.gz --wildcards "*.fasta.gz" --exclude "*_DNA*"

cat Bacteria/*/*.fasta.gz Archaea/*/*.fasta.gz \
    Eukaryota/*/*.fasta.gz Viruses/*/*.fasta.gz \
    > reference_proteomes.fasta.gz
