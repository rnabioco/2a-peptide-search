
#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO

import csv
import sys
import gzip
import pdb

for row in csv.DictReader(gzip.open("all-sequences.tsv.gz", "rt"), delimiter="\t"):
    pdb.set_trace()

for record in SeqIO.parse(gzip.open(sys.argv[1], "rt"), "fasta"):
    pdb.set_trace()

for aln in AlignIO.parse(sys.argv[2], "stockholm"):
    for record in aln.alignment.sequences:
        pass

