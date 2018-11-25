#!/usr/bin/env python
from Bio import SeqIO

def split_fasta(infile):
    # Split fasta file into 10 smaller filesfiles for each sequence
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

