#!/usr/bin/python
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict
import tqdm

# To Run script, specify:
# targets_filename - fasta file with targets sequences
# antismash_filename - output file of ks_target tandem search
# gbdir1 and gbdir2 - folders with genbank files
# outdir - outdir where output will be written

target_to_seq = {}

# Make a dictionary for target ids and their sequences
targets_filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/targets.616.fa"

for record in SeqIO.parse(open(targets_filename, "rU"), "fasta"):
    targetid = record.id
    targetseq = record.seq
    target_to_seq[targetid] = targetseq
    print targetid, targetseq[:10]

# Read antismash output KS fasta file
# ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
antismash_ksfasta = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/KS.616.10kb.fasta.cdhit.90"
gbid_to_target = defaultdict(list)

for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
    gbidfull = record.id
    gbidfull = gbidfull.strip()
    gbid, targetid = gbidfull.split("|")[:2]
    print gbid, targetid
    gbid_to_target[gbid].append((targetid))

for gbid in tqdm.tqdm(sorted(gbid_to_target.keys())):
    targetids = gbid_to_target[gbid]
    for targetid in targetids:
        # if targetid == "mupM_Ile-tRNA-syn" and gbid == "CP025542":
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

        outdir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Second_copy/blast_out"
        blast_evalue_cutoff = 1e-8
        outfilename = os.path.join(outdir, targetid + "." +
                                   str(blast_evalue_cutoff) +
                                   "." + gbid + ".out")
        # print outfilename
        print queryfile
        print gbfile
        # tblastn -query targets_fasta/mupM_Ile-tRNA-syn -subject /mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/assembly_gb/CP025542.gbff -out out.txt
        subprocess.call(["tblastn", "-query", queryfile,
                         "-subject", gbfile,
                         "-out", outfilename,
                         "-evalue", str(blast_evalue_cutoff),
                         "-outfmt",  "6 qseqid sseqid sstart send nident qlen slen evalue"])
