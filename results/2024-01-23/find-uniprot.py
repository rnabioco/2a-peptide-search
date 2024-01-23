#! /usr/bin/env python

import sys
import gzip
from collections import defaultdict

lookup = defaultdict(set)

for line in gzip.open(sys.argv[1], 'rt'):

    if line.startswith("#"): continue
    fs = [f.strip() for f in line.split("\t")]
    
    name, seq = fs[0], fs[1] 
    lookup[seq].add(name)

# check for unique
for seq, names in lookup.items():
    if len(names) > 1:
        print(len(names), names, sep = "\t")


