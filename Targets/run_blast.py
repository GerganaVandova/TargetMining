#!/usr/bin/env python
# Blast putative target query genes against nucleotide biosynthetic clusters databases

import os
import subprocess
import glob
import sys
import math

# Blast parameters
# num_alignments = 100000
blast_evalue_cutoff = 1e-50
# GV not used in this script; added in parse_script
# blast_evalue_cutoff_parse = 1
blastoutdir = "blast_results"
subprocess.call(["mkdir", "-p", blastoutdir])

# There is a database for each cluster
blastdbdir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Targets/blastdb/"
mibigids = []
f = open('mibig_ids.txt', "r")
mibig_ids = f.readlines()
for mibig_id in mibig_ids:
    mibig_id = mibig_id.strip()
    mibigids.append(mibig_id)
    # print mibig_id

# Specify which queries to blast (KS*.seq, CLF*.seq, etc)
blast_query_dir = "uniprot_fasta/"
blast_query_files = glob.glob("uniprot_fasta/*.fasta")
print blast_query_files
blast_query_names = map(os.path.basename, blast_query_files)
print blast_query_names

# Blast KS, etc domains:
for blast_query_name in blast_query_names:
    for mibigid in mibigids:
        dbdir = os.path.join(blastdbdir, mibigid)
        print dbdir
        blast_query_file = os.path.join(blast_query_dir, blast_query_name)
        print blast_query_name, mibigid, "\n"
        outfilename = os.path.join(blastoutdir, blast_query_name + "." + str(blast_evalue_cutoff) + "." + mibigid + ".out")
        
        #if os.path.exists(outfilename) == True:
        #    print outfilename, " exist"
        #    continue

        outfilename = os.path.join(blastoutdir, blast_query_name + "." + str(blast_evalue_cutoff) + "." + mibigid + ".out")
        print outfilename
        subprocess.call(["blastp", "-query", blast_query_file, "-db", dbdir, "-out",  outfilename, "-evalue", str(blast_evalue_cutoff), "-outfmt", "6"])
