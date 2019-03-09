#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
# import tqdm

# Get phyla from taxa file
# BDBI01000023   Nocardia sp.
gbid_to_species = {}
taxafilename = "../Genbank/species.txt"
taxafile = open(taxafilename).readlines()
for line in taxafile:
    gbid, name = line.split("\t")
    gbid_to_species[gbid] = name

# Get phyla from taxa file
# BDBI01000023   Bacteria    Actinobacteria  Corynebacteriales   Nocardiaceae    Nocardia
gbid_to_phyla = {}
taxafilename = "../Genbank/taxa.txt"
taxafile = open(taxafilename).readlines()
for line in taxafile:
    feats = line.split("\t", 3)
    if len(feats) < 3:
        continue
    gbid, taxa, phyla, rest = feats
    gbid_to_phyla[gbid] = phyla

f = open("KS.609.5kb.fasta.descr", "w")
# f1 = open("KS.609.5kb.fasta", "w")

# Read KS fasta file and write phyla for each gbid
# >KB907307|DEG10180214_Oligopeptide_transport_ATP-binding_protein_oppD|1954686|1955110|1951543|1953327|cluster-1|terpene-t1pks|1918076-2020118|1359
for record in SeqIO.parse(open("KS.609.5kb.fasta.origlongnames", "rU"), "fasta"):
    gbidfull = record.id
    seq = record.seq
    gbid, target, _, _, _, _, _, cl_type = gbidfull.split("|")[:8]
    d1 = gbidfull.split("|")[2:]
    d = "|".join(d1)
    target_short, target_descr = target.split("_", 1)
    gbid_new = "|".join([gbid, target_short, d])
    print "d\t%s\ntarget_short\t%s\ntarget_descr\t%s\ngbid_new\t%s" % \
          (d, target_short, target_descr, gbid_new)

    # f1.write(">%s\n%s\n" % (gbid_new, str(seq)))

    if gbid not in gbid_to_phyla.keys():
        gbid_to_phyla[gbid] = "None"
    if gbid not in gbid_to_species.keys():
        gbid_to_species[gbid] = "None"

    descr1 = "_".join([gbid, target, cl_type,
                       gbid_to_phyla[gbid],
                       gbid_to_species[gbid]])

    f.write("%s\t%s\n" % (gbid_new, descr1))
    print "%s\t%s" % (gbid_new, descr1)

f.close()
# f1.close()
