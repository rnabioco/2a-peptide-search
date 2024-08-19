#! /usr/bin/env bash

awk '{print ">"$1"\n"$2}' < seqs.tsv > seqs.faa

hmmsearch --max --tblout seqs.class-1.tab ../../curated-models/2A-class-1.hmm seqs.faa > seqs.class-1.hmmsearch
hmmsearch --max --tblout seqs.class-2.tab ../../curated-models/2A-class-2.hmm seqs.faa > seqs.class-2.hmmsearch
