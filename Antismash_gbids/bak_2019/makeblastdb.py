#!/usr/bin/env python
from Bio import SeqIO
import os

seq_dir = "seq_dir/"
for seqfile in os.listdir(seq_dir):
    f = os.path.join(seq_dir, seqfile)
    db = "antismashdb_" + seqfile
    makeblastdb = "makeblastdb -in " + f + " -dbtype prot -out " + db
    print makeblastdb
    os.system(makeblastdb)

