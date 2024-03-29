#!/usr/bin/env python
# Blast query genes against nucleotide databases

import os
import subprocess
import glob

blast_query_dir = "Manual"

# updated blast databases on 01/27/2016
blastdbdir = "/mnt/gnpn/gnpn/projects/orphanpks/blast/db"

# Blast parameters
num_alignments = 100000
blast_evalue_cutoff = 1
# GV not used in this script; added in parse_script
blast_evalue_cutoff_parse = 1
blastoutdir = "blast_results"


subprocess.call(["mkdir", "-p", blastoutdir])

dbs = ["nt", "refseq_genomic", "other_genomic", "env_nt", "patnt", "sts", "htgs", "tsa_nt", "wgs", "gss", "est_others", "asembliesntdb"]
#dbs = ["asembliesntdb"]

#Specify which queries to blast (KS*.seq, CLF*.seq, etc)
blast_query_files = glob.glob("%s/KS*.seq" % blast_query_dir)
blast_query_names = map(os.path.basename, blast_query_files)

# Blast KS, etc domains:
for db in dbs:
    for blast_query_name in blast_query_names:
        dbdir = os.path.join(blastdbdir, db)
        blast_query_file = os.path.join(blast_query_dir, blast_query_name)
        print blast_query_name, db, "\n"
        outfilename = os.path.join(blastoutdir, blast_query_name + "." + str(num_alignments) + ".alignments.evalue" + str(blast_evalue_cutoff) + "." + db)
        print outfilename
        blastcmd = "tblastn", "-query", blast_query_file, "-db", dbdir, "-out",  outfilename, "-num_alignments", str(num_alignments), "-seg",  "no",  "-evalue", str(blast_evalue_cutoff),  "-num_threads", "64"
        print blastcmd
        subprocess.call(["tblastn", "-query", blast_query_file, "-db", dbdir, "-out",  outfilename, "-num_alignments", str(num_alignments), "-seg",  "no",  "-evalue", str(blast_evalue_cutoff),  "-num_threads", "64"])


# parse_blast output:
#	print OUT "$hitname\t$start\t$end\t$hitdesc\t$frame\t$strand\t$hitlen\t$evalue\t$hitseq\n";

parse_script = "scripts/getGBHitCoordsFilesInDir.pl"
parse_blast_dir = "blast_results_parse"

subprocess.call(["mkdir", "-p", parse_blast_dir])
subprocess.call(["perl", parse_script, blastoutdir, parse_blast_dir])
