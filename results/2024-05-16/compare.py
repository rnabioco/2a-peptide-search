
#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
from Bio import SearchIO

from collections import defaultdict

import csv
import sys
import gzip
import pdb

ref_seqs = dict()
row_num = 1
for row in csv.DictReader(gzip.open("all-sequences.tsv.gz", "rt"), delimiter="\t"):
    # id in the tsv, hmmer adds underscores
    name = row['Name'].replace(" ", "_")
    anal_id = row['Analysis ID']
    # id in the fasta
    name_fa = f"{name}|{row_num}|{anal_id}".replace('N/A', 'NA')
    

    seq = row['Sequence']
    status = row['Old/new']

    ref_seqs[name_fa] = (name, seq, status)

    row_num += 1

hit_records = defaultdict(set)
for result in SearchIO.parse("lit-annotations.class-1.hmmsearch.tab", "hmmer3-tab"):
    for hit in result.hits:
        if hit.domain_included_num == 1:
            hit_records["class-1"].add(hit.id)
for result in SearchIO.parse("lit-annotations.class-2.hmmsearch.tab", "hmmer3-tab"):
    for hit in result.hits:
        if hit.domain_included_num == 1:
            hit_records["class-2"].add(hit.id)

for ref_name, ref_info in ref_seqs.items():
    for class_name, hits in hit_records.items():
        if ref_name in hits:
            print(ref_name, *ref_info, class_name, "found", sep = "\t")
        else:
            print(ref_name, *ref_info, class_name, "absent", sep = "\t")

