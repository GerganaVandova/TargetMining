#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
import tqdm

gbid_to_phyla = {}

# Get phyla from taxa file
#BDBI01000023   Bacteria    Actinobacteria  Corynebacteriales   Nocardiaceae    Nocardia
taxafilename = "../Genbank/taxa.txt"
taxafile = open(taxafilename).readlines()
for line in taxafile:
    feats = line.split("\t")
    if len(feats) >= 3:
        gbid, taxa, phyla = feats[:3]
        print gbid, taxa, phyla
        gbid_to_phyla[gbid] = phyla


#f = open("out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.fasta.phyla", "w")
f = open("out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.withFabF.fasta.cdhit.90.phyla", "w")


# Read KS fasta file and write phyla for each gbid
#ACXX02000001_38012-39229
#fasta_file = "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.fasta"
fasta_file = "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.withFabF.fasta.cdhit.90"
for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
    gbidfull = record.id
    if "FabF" in gbidfull:
        continue
    #print gbidfull
    gbid, coord  = gbidfull.rsplit("_", 1)
    if gbid in gbid_to_phyla.keys():
        print gbidfull, gbid_to_phyla[gbid]
        f.write("%s\t%s\n" % (gbidfull, gbid_to_phyla[gbid]))
    else:
        print gbidfull, None
        f.write("%s\tNone\n" % (gbidfull))

f.close()

