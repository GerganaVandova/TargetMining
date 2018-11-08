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


def split_mibig_fasta(infile):
    # Split mibig fasta file into subfiles for each cluster
    with open(infile, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        for i in xrange(len(records)):
            print(records[i].id)
            mibigid = records[i].id.split('|')[0]
            print mibigid
            outfilename = 'mibig_fasta/' + mibigid + '.fasta'
            f = open(outfilename, 'wa')
            f.write(">%s\n" % mibigid)
            f.write(str(records[i].seq))
            f.close()


def split_fasta(infile):
    # Split fasta file into files for each sequence
    # Store seq files in uniprot_seqs/ folder

    with open(infile, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        for record in records:
            new_name = record.name.replace('|', '_')
            outfilename = 'uniprot_seqs/' + new_name + '.seq'
            print outfilename
            f = open(outfilename, 'wa')
            f.write(">%s\n" % new_name)
            f.write(str(record.seq))
            f.close()



# infile = 'MIBiG_prot_seqs_1.4.fasta'
infile = 'all_uniprot.fa'
split_fasta(infile)
