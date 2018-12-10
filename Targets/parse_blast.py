#!/usr/bin/env python
# Blast putative target query genes against nucleotide biosynthetic clusters databases

import os
import subprocess
import glob
import sys
import math


# BGC0000185  tartrolon   Polyketide  None    
id_to_name = {}

target_to_cluster = {}


f = open("mibig_clusters.txt", 'r')
for line in f.readlines():
    line = line.strip()
    params = line.split("\t")
    mibigid = params[0]
    name = params[1]
    id_to_name[mibigid] = (params[1:])


blast_files = glob.glob("blast_results/*.out")

for blast_file in blast_files:
    q = blast_file.split(".")[0]
    target = q.split('blast_results/')[1]
    f = open(blast_file, 'r')
    for line in f.readlines():
        line = line.strip()
        
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")
        pident = float(pident)
        mibigid = sseqid.split("|")[0]
    
        if pident > 50:
            if mibigid == "BGC0001482":
                print "%s\t%s\tNoparams\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mibigid, target, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
            else:
                target_to_cluster[(mibigid, target)] = (sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

sorted(target_to_cluster.keys())

for key in target_to_cluster.keys():
    mibigid, target = key
    sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = target_to_cluster[key]
    if len(id_to_name[mibigid]) == 4:
        a, b, c, d = id_to_name[mibigid]
#        if "Antibacterial" in d or "Cytotoxic" in d or c != None or c != "Unknown":
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mibigid, target, a, b, c, d, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

    elif len(id_to_name[mibigid]) == 3:
        a, b, c = id_to_name[mibigid]
        print "%s\t%s\t%s\t%s\t%s\tNoparams\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mibigid, target, a, b, c, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

    else:
        a, b = id_to_name[mibigid]
        print "%s\t%s\t%s\t%s\tNoparams\tNoparams\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mibigid, target, a, b, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
#print len(target_to_cluster)


