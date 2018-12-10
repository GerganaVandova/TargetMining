#!/usr/bin/env python
# Blast putative target query genes against nucleotide biosynthetic clusters databases

import os
import subprocess
import glob
import sys
import math

f = open("all.out", 'r')
for line in f.readlines():
    line = line.strip()

    #query_name ca2_atpase
    #qseqid sp|A0R1E8|PKS5_MYCS2
    #sseuqid BGC0001795|c1|46569-59957|+|no_locus_tag|malonyl_CoA-acyl_carrier_protein_transacylase|ATX68126.1
    #mibigid BGC0001795

    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split("\t")
    pident = float(pident)
    mibigid = sseqid.split("|")[0]
    print mibigid, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore

