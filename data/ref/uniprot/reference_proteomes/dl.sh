#! /usr/bin/env bash

#BSUB -J dload-ref-proteomes
#BSUB -eo logs/%J.errout

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2024_02.tar.gz
