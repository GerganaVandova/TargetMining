#!/usr/bin/python
import sys
import json
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict

target_to_seq = {}

# Make a dictionary for target ids and their sequences
targets_filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/targets.12.fa"

for record in SeqIO.parse(open(targets_filename, "rU"), "fasta"):
    targetid = record.id
    targetseq = record.seq
    # print targetid, targetseq[:10]
    target_to_seq[targetid] = targetseq


# Read antismash output file
# ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
antismash_filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/out.12.filtered.10kb"
antismash_file = open(antismash_filename).readlines()

gbid_to_target = defaultdict(list)

for line in antismash_file:
    line = line.strip()
    features = line.split("|")
    gbid, targetid = features[:2]
    print gbid, targetid
    gbid_to_target[gbid].append((targetid))

for gbid in sorted(gbid_to_target.keys()):
    targetids = gbid_to_target[gbid]
    for targetid in targetids:
        targetseq = target_to_seq[targetid]
        print gbid, targetid, targetseq[:10]

        gbdir1 = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbdir"
        gbdir2 = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/assembly_gb"

        path1 = os.path.join(gbdir1, gbid + ".gb")
        path2 = os.path.join(gbdir2, gbid + ".gbff")

        if os.path.exists(path1) is True:
            gbfile = path1
        elif os.path.exists(path2) is True:
            gbfile = path2
        else:
            print "filename not correclty parsed"
            print gbfile1
            print gbfile2

        queryfile = "targets_fasta/%s.fasta" % targetid
        # print queryfile
        # print gbfile

        outdir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Second_copy"
        blast_evalue_cutoff = 1e-8
        outfilename = os.path.join(outdir, targetid + "." +
                                   str(blast_evalue_cutoff) +
                                   "." + gbid + ".out")
        # print outfilename

        # EF-Tu.1e-08.NZ_KQ948231.out
        # if targetid == "borI_Thr-tRNA-syn" and gbid == "KT362046":
        #     print queryfile
        #     print gbfile
        # tblastn -query targets_fasta/EF-Tu.fasta -subject /mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/assembly_gb/KB899005.gbff -out out.txt
        subprocess.call(["tblastn", "-query", queryfile,
                         "-subject", gbfile,
                         "-out", outfilename,
                         "-evalue", str(blast_evalue_cutoff),
                         "-outfmt",  "6 qseqid sseqid sstart send nident qlen slen evalue"])
