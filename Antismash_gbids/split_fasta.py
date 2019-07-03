#!/usr/bin/env python
from Bio import SeqIO
import json
import os
from pprint import pprint
import glob
import codecs
import sys
import requests
import urllib
import subprocess
UTF8Writer = codecs.getwriter('utf8')
sys.stdout = UTF8Writer(sys.stdout)


def get_targets(antismash_ksfasta):
    targets = set()
    for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
        gbidfull = record.id
        target = gbidfull.split("|")[1]
        targets.add(target)
    return(targets)


# For each target make a fasta file with all KS sequences
def split_fasta(infile, targets):

    for target in targets:
        outfilename = "fasta/%s.fasta" % target
        f = open(outfilename, 'a')
        with open(infile, "rU") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            for record in records:
                name = record.name
                qtarget = name.split("|")[1]
                if qtarget == target:
                    f.write(">%s\n" % name)
                    f.write(str(record.seq))
                    f.write("\n")

fasta_file = 'KS.14.10kb.fasta'
targets = get_targets(fasta_file)
# split_fasta(fasta_file, targets)

# # Run cdhit:
# for fasta_file in os.listdir("fasta/"):
#     out_file = fasta_file + ".cdhit.90"
#     print fasta_file, out_file
#     # run_cdhit(dir, fasta_file, outfile)
#     cmd = "~maureenh/Cdhit/cd-hit-v4.6.1-2012-08-27/cd-hit -i " + "fasta/" + fasta_file + " -o " + "fasta/" + out_file + " -c .9 -d 200 -M 10000"
#     print cmd
#     os.system(cmd)

for target in targets:
    cdhitfile = "KS.14.10kb.fasta.cdhit.90t"
    cmd = "cat KS.14.10kb.fasta.cdhit.90t | grep " + target + " |wc"
    print cmd
    os.system(cmd)
