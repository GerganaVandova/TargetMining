#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict

# Get target name from targets fasta file
# >DEG10180001_Molybdopterin_biosynthesis_mog_protein
gbid_to_names = {}
targetfilename = "../Antismash_gbids/targets.616.fa.cleannames"
for record in SeqIO.parse(open(targetfilename, "rU"), "fasta"):
    gbidfull = record.id
    gbid, name = gbidfull.split("_", 1)
    name = name.replace("/", "_")
    name = name.replace(" ", "_")
    name = name.replace("\n", "_")
    gbid_to_names[gbid] = name

# Get species name from species file
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

f = open("KS.616.10kb.fasta.filtered.descr", "w")
# Read KS fasta file and write phyla for each gbid
# KB907307|DEG10180214|1954686|1955110|1951543|1953327|cluster-1|terpene-t1pks|1918076-2020118|1359
for record in SeqIO.parse(open("KS.616.10kb.fasta.cdhit.90", "rU"), "fasta"):
    gbidfull = record.id
    seq = record.seq
    gbid, target_id, rest = gbidfull.split("|", 2)
    cl_num = gbidfull.split("|")[6]
    cl_type = gbidfull.split("|")[7]
    # print gbid, target_id, rest
    if gbid not in gbid_to_phyla.keys():
        gbid_to_phyla[gbid] = "None"
    if gbid not in gbid_to_species.keys():
        gbid_to_species[gbid] = "None"

    # With phyla
    # descr = "_".join([gbid, target_id, cl_type,
    #                   cl_num,
    #                   gbid_to_phyla[gbid],
    #                   gbid_to_names[target_id],
    #                   gbid_to_species[gbid]])

    # Without phyla nd target name
    descr = "_".join([gbid, target_id, cl_type,
                    cl_num,
                    gbid_to_names[target_id],
                    gbid_to_species[gbid]])

    f.write("%s\t  %s\n" % (gbidfull, descr))
    print "%s\t  %s" % (gbidfull, descr)

f.close()
