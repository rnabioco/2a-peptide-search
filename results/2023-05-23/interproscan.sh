#! /usr/bin/env bash

#BSUB -J interproscan[1-2]
#BSUB -e logs/interproscan.%J.%I.err
#BSUB -o logs/interproscan.%J.%I.out
#BSUB -n 16

interproscan=$HOME/src/interproscan-5.62-94.0/interproscan.sh
i=$LSB_JOBINDEX
project=$HOME/devel/rnabioco/2a-peptide-search
fa=$project/results/2023-04-22/class-$i.fa

$interproscan \
    -appl FunFam,SFLD,PANTHER,Gene3D,Hamap,PRINTS,Coils,SUPERFAMILY,SMART,CDD,PIRSR,AntiFam,Pfam,MobiDBLite,PIRSF,NCBIfam \
    -i $fa \
    -f tsv \
    -dp \
    -o interproscan.class-$i.tsv \
    --cpu 16

