#!/usr/bin/env python
from Bio import SeqIO
import os

blastdb_dir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb"
blastdb_names = ["/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.89k.coord.fasta"]

for blastdb_name in blastdb_names:
    outfile1 = "out.targets.9" + blastdb_name.split("89k.coord")[1]
    outfile = outfile1.split(".fasta")[0]
    blastp = "blastp -db " + blastdb_name + " -query targets.9.fa -outfmt " + '"6 qseqid sseqid sstart send nident qlen slen evalue"' +  " -out " + outfile
    print blastp
    os.system(blastp)
