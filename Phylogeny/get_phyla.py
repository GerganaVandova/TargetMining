#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
import tqdm

gbid_to_phyla = {}

# Get phyla from taxa file
# BDBI01000023   Bacteria    Actinobacteria  Corynebacteriales   Nocardiaceae    Nocardia
taxafilename = "../Genbank/taxa.txt"
taxafile = open(taxafilename).readlines()
for line in taxafile:
    feats = line.split("\t", 3)
    if len(feats) < 3:
        continue
    gbid, taxa, phyla, rest = feats
    gbid_to_phyla[gbid] = phyla

f = open("KS.12.20kb.fasta.phyla", "w")

# Read KS fasta file and write phyla for each gbid

# >CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
for record in SeqIO.parse(open("KS.12.20kb.fasta", "rU"), "fasta"):
    gbidfull = record.id
    gbid, rest = gbidfull.split("|", 1)
    coord = gbidfull.rsplit("|", 2)[1]
    if gbid not in gbid_to_phyla.keys():
        f.write("%s\tNone\n" % gbidfull)
        print "%s\tNone" % gbidfull
        continue
    f.write("%s\t%s\n" % (gbidfull, gbid_to_phyla[gbid]))
    print "%s\t%s" % (gbidfull, gbid_to_phyla[gbid])

f.close()

# KS.12.10kb.fasta.phyla
# PYAX01000012|120458-199012|Actinobacteria
