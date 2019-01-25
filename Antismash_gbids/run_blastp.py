#!/usr/bin/env python
from Bio import SeqIO
import os
import ntpath

targets_file = "targets.12.fa"
dbdir = "antismashdb"
blastdb_names = [os.path.join(dbdir, "cluster_genes.89k"),
                 os.path.join(dbdir, "cluster_genes.21k")]

print blastdb_names
for blastdb_name in blastdb_names:
    outfile = "out." + targets_file + "." + ntpath.basename(blastdb_name)
    print outfile
    blastp = "blastp -db " + blastdb_name + \
             " -query " + targets_file + \
             " -outfmt " + '"6 qseqid sseqid sstart send nident qlen slen evalue"' + \
             " -out " + outfile
    print blastp
    os.system(blastp)
