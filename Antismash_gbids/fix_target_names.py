#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
f = open("targets.609.fa.fixed", "w")
# f1 = open("KS.609.5kb.fasta", "w")


# Read KS fasta file and write phyla for each gbid
# >DEG10180001_Molybdopterin_biosynthesis_mog_protein
for record in SeqIO.parse(open("targets.609.fa", "rU"), "fasta"):
    gbidfull = record.id
    seq = record.seq
    gbid, rest = gbidfull.split("_", 1)
    print gbid
    f.write(">%s\n%s\n" % (gbid, seq))

f.close()
