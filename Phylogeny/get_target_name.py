#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict

f = open("KS.609.5kb.fasta.target", "w")

# Read KS fasta file and write phyla for each gbid
# >CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
for record in SeqIO.parse(open("KS.609.5kb.fasta", "rU"), "fasta"):
    gbidfull = record.id
    target = gbidfull.split("|")[1]
    f.write("%s\t%s\n" % (gbidfull, target))
    print "%s\t%s" % (gbidfull, target)

f.close()
