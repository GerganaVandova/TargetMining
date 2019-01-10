#!/usr/bin/python
import sys
import json
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict


#./filter_second_copy.py > blast_out.second_copy.filtered.pident.30.genome_len
#    outfile: blast_out.second_copy.filtered.pident.30
#             blast_out.second_copy.filtered.pident.30.genome_len


PIDENTCUTOFF = 30
PIDENTCUTOFF_FAB = 60

outfiltered_file = open("blast_out.second_copy.filtered.pident.30", "w")
blastdir = "blast_out/"

pairs = defaultdict(int)

# AdmT_ACC.1e-08.KB894406.out blast_out/AdmT_ACC
# AdmT_ACC  Subject_1   48.02   205138  205893  121 304 295743  4e-70 
for blastoutfile in os.listdir(blastdir):
    f = os.path.join(blastdir, blastoutfile)
    target, _, cluster = f.split(".")[:3]
    #print blastoutfile, target, cluster
    lines = open(f).readlines()
    for line in lines:
        target, _, pident, sstart, send, nident, qlen, slen, evalue = line.strip().split("\t")
        cutoff = PIDENTCUTOFF
        if target == "PtmP3_FabB-F":
            cutoff = PIDENTCUTOFF_FAB
        if float(pident) < cutoff:
            continue
        s = "\t".join([target, cluster, pident, sstart, send, nident, qlen, slen, evalue])
        outfiltered_file.write("%s\n" %s)
        pairs[(target, cluster)] += 1


# Read genome_lengths file 
# A70054    False   1223 
gbid_to_len = {}
genome_lengths_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/genome_lengths.txt"
lines = open(genome_lengths_file, "r").readlines()
for line in lines:
    gbid, complete_genome, length = line.strip().split("\t")
    gbid_to_len[gbid] = ((complete_genome, length))


# Read antismash output file 
#AdmT_ACC   ACXX02000001    1   14512-116691    transatpks-nrps Cpap_3701_acetyl-CoA    31720   32586   118.0   304.0   288 0.39    1e-67   38012   39229   5426    102179
antismash_filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"
antismash_file = open(antismash_filename).readlines()[1:]
for line in antismash_file:
    line = line.strip()
    features = line.split("\t")
    targetid, gbid  = features[:2]
    complete_genome, length = gbid_to_len[gbid]
    length = int(length)
    if complete_genome == "False":
        if length > 4000000:
            complete_genome = "probably true"
    print "%s\t%s\t%s\t%s\t%s" % (targetid, gbid, pairs[(targetid, gbid)], length, complete_genome)

