#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict

# Get target name from targets fasta file
# >DEG10180001_Molybdopterin_biosynthesis_mog_protein
gbid_to_names = {}
targetfilename = "../Antismash_gbids/targets.609.fa.longnames"
for record in SeqIO.parse(open(targetfilename, "rU"), "fasta"):
    gbidfull = record.id
    gbid, name = gbidfull.split("_", 1)
    name = name.replace("/", "_")
    gbid_to_names[gbid] = name

f = open("KS.609.5kb.fasta.filtered", "w")

for record in SeqIO.parse(open("KS.609.5kb.fasta", "rU"), "fasta"):
    gbidfull = record.id
    seq = record.seq
    gbid, target_id, rest = gbidfull.split("|", 2)
    cl_type = gbidfull.split("|")[6]
    target_name = gbid_to_names[target_id]
    # print gbid, target_id, target_name

    if "[acyl-carrier-protein]" in target_name:
        print "skipping %s" % target_name
        continue
    if "transport" in target_name:
        print "skipping %s" % target_name
        continue

    f.write(">%s\n%s\n" % (gbidfull, str(seq)))
    # print "%s\n%s\n" % (gbidfull, seq[:10])

f.close()
