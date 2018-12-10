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


# For each cluster make a fasta file with all proteins 

def split_fasta(infile):
    # Split fasta file into files for each sequence
    # Store seq files in uniprot_seqs/ folder

    with open(infile, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        for record in records:
            name = record.name
            #>BGC0000001|c1|1-1083|+|no_locus_tag|protein_methyltransferase|AEK75490.1
            print name
            cluster = name.split("|")[0]
            outfilename = "clusters_fasta/%s.fasta" % cluster
            f = open(outfilename, 'a')
            f.write(">%s\n" % name)
            f.write(str(record.seq))
            f.write("\n")

infile = 'MIBiG_prot_seqs_1.4.fa'
#infile = 'all_uniprot.fa'
split_fasta(infile)
