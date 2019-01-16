#!/usr/bin/env python
from Bio import SeqIO
import os
import sys

seq_dirs = ["ks_seqs/", "target_seqs/"]
for seq_dir in seq_dirs:
    for seqfile in os.listdir(seq_dir):
        print seqfile
        f = os.path.join(seq_dir, seqfile)
        dbname = seqfile.split(".fasta")[0]
        makeblastdb = "makeblastdb -in " + f + " -dbtype prot -out " + dbname
        print makeblastdb
        os.system(makeblastdb)
