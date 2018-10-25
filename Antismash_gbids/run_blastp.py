#!/usr/bin/env python
from Bio import SeqIO
import os

blastdb_dir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb"
blastdb_names = ["/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.00.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.01.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.02.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.03.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.04.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.05.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.06.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.07.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.08.fasta",
                 "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismashdb/antismashdb_sequences.faa.95k.coord.09.fasta"]

for blastdb_name in blastdb_names:
    outfile1 = "out.targets.609" + blastdb_name.split("84k.coord")[1]
    outfile = outfile1.split(".fasta")[0]
    blastp = "blastp -db " + blastdb_name + " -query targets.609.fa -outfmt " + '"6 qseqid sseqid sstart send nident qlen slen evalue"' +  " -out " + outfile
    print blastp
    os.system(blastp)
