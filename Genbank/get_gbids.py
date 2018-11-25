#!/usr/bin/env python
from Bio import SeqIO

filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
nt_outfile = "gbids.nt.txt"
assembly_outfile = "gbids.assembly.txt"

nt_f = open(nt_outfile, 'a')
assembly_f = open(assembly_outfile, 'a')

"""
for record in SeqIO.parse(filename, "fasta"):
    recordid = record.id
    if "__" in recordid:
        gbid = recordid.split("__")[0]
        assembly_f.write(gbid)
        assembly_f.write("\n")
    else:
        gbid = recordid.split('.')[0]
        nt_f.write(gbid)
        nt_f.write("\n")
    print gbid
"""

gbids = []
nt_gbids = []
assembly_gbids = []

for record in SeqIO.parse(filename, "fasta"):
    recordid = record.id
    gbid = recordid.split('.')[0]
    gbids.append(gbid)

for gbid in gbids:
    if "__" in gbid:
        new_gbid = gbid.split("__")[0]
        assembly_gbids.append(new_gbid)
        assembly_f.write(new_gbid)
        assembly_f.write("\n")
        print new_gbid
    else:
        nt_gbids.append(gbid)
        nt_f.write(gbid)
        nt_f.write("\n")
        print gbid
