#! /usr/bin/env bash

#BSUB -J extract-ref

tar zxvf Reference_Proteomes_2024_02.tar.gz --wildcards "*.fasta.gz" --exclude "*_DNA*"
