#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict

gbid_to_species = {}

# Get phyla from taxa file
# BDBI01000023   Nocardia sp.
taxafilename = "../Genbank/species.txt"
taxafile = open(taxafilename).readlines()
for line in taxafile:
    gbid, name = line.split("\t")
    gbid_to_species[gbid] = name

f = open("KS.12.20kb.fasta.descr.species", "w")

# Read KS fasta file and write phyla for each gbid
# >CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
for record in SeqIO.parse(open("KS.12.20kb.fasta", "rU"), "fasta"):
    gbidfull = record.id
    gbid, target, _, _, _, _, _, cl_type = gbidfull.split("|")[:8]
    # if gbid not in gbid_to_species.keys():
    #     f.write("%s\tNone\n" % gbidfull)
    #     print "%s\tNone" % gbidfull
    #     continue
    descr = "_".join([gbid, target, cl_type, gbid_to_species[gbid]])
    f.write("%s\t%s\n" % (gbidfull, descr))
    print "%s\t%s" % (gbidfull, descr)

f.close()
