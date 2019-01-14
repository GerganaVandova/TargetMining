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


def split_fasta(infile):
    # Split fasta file into files for each sequence
    # Store seq files in uniprot_seqs/ folder

    with open(infile, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        for record in records:
            name = record.name
            outfilename = "fasta/fasta_split/%s.fasta" % name
            f = open(outfilename, 'a')
            f.write(">%s\n" % name)
            f.write(str(record.seq))
            f.write("\n")

infiles = ['fasta/KS_seq/all.KS.fasta', 'fasta/target_seq/all.target.fasta']

for infile in infiles:
    split_fasta(infile)
