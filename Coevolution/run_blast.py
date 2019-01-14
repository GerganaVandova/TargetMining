#!/usr/bin/python
import sys
import json
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict
import os.path


def get_targets(antismash_outfilename):
    # Get list of target gene names

    targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target = line.split("\t")[0]
        targets.add(target)
    return(targets)


def paiwise_blast(blastdb_dir, fasta_dir, target_name):
    # For each target, do paiwise blasts for KSs and targets
    # from the same clusters
    for fasta_seq in os.listdir(fasta_dir):
        fasta_seq_basename = fasta_seq.split(".fasta")[0]
        outfile = target_name + ".out"
        if fasta_seq_basename != target_name:
            continue

        for blastdb_filename in os.listdir(blastdb_dir):
            # print blastdb_filename
            blastdb = blastdb_filename.rsplit(".", 1)[0]
            # print blastdb
            if blastdb != target_name:
                continue

            if os.path.isfile(outfile):
                continue

            blastp = "blastp -db " + os.path.join(blastdb_dir, blastdb) \
                + " -query " + os.path.join(fasta_dir, fasta_seq) \
                + " -outfmt " \
                + '"6 qseqid sseqid sstart send nident qlen slen evalue"' \
                + " -out " + outfile
            # blastp -db blastdb/GyrB-R -query GyrB-R.KS.fasta -outfmt "6 qseqid sseqid sstart send nident qlen slen evalue" -evalue 1e-8 -out GyrB-R.KS.out
            print blastp

            os.system(blastp)


blastdb_dir = "blastdb"
fasta_dir = "fasta_seqs"
antismash_outfilename = \
    "../Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"


targets = get_targets(antismash_outfilename)
targets = ["GyrB-R"]
for target in targets:
    target_ks = target + ".KS"
    paiwise_blast(blastdb_dir, fasta_dir, target)
    paiwise_blast(blastdb_dir, fasta_dir, target_ks)
