#! /usr/bin/env python

import sys
import gzip

from pyfaidx import Fasta
from collections import defaultdict

lookup = defaultdict(set)

for line in gzip.open(sys.argv[1], 'rt'):

    if line.startswith("#"): continue
    fs = [f.strip() for f in line.split("\t")]
    
    name, seq, status = fs[0], fs[1], fs[3]

    if status == 'Old':
        lookup[seq].add(name)

# check for unique
#for seq, names in lookup.items():
#    if len(names) > 1:
#        print(len(names), names, sep = "\t")

fasta = Fasta(sys.argv[2])
fasta_ids = fasta.keys()
for fasta_id in fasta_ids:
    for seq in lookup:
        if seq in fasta[fasta_id]:
            print(fasta_id, seq, sep = "\t")

