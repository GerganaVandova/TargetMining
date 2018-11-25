#!/usr/bin/env python
import subprocess
import os
import os.path
import tqdm
from Bio import Entrez
from multiprocessing import Pool

TEMP_ANTISMASH_OUTPUT = "/home/gvandova/antismash_output_2018/"
FINAL_ANTISMASH_OUTPUT = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output"
GB_FOLDER = "/home/gvandova/gbdir"


def f(gbid):
    #print "\n", "***********", "\n", gbid
    temp_dir = os.path.join(TEMP_ANTISMASH_OUTPUT, gbid)
    final_dir = os.path.join(FINAL_ANTISMASH_OUTPUT, gbid)

    local_filename = os.path.join(GB_FOLDER, gbid + ".gb")
    geneclusters_filename = os.path.join(final_dir, "geneclusters.js")
    print geneclusters_filename
    if os.path.exists(geneclusters_filename) == True:
        print geneclusters_filename, " exist"
        return
    if os.path.exists(final_dir):
        print final_dir, "exists"
	return
    #antismash = "python /home/gvandova/antismash/run_antismash.py" + " --outputfolder " + \
    #            antismashoutdir + " " + local_filename
    antismash = "/home/gvandova/bin/run_antismash" + " " + local_filename + " " + temp_dir
    print antismash
    print "start %s" % gbid
    os.system(antismash)
    subprocess.call(["rm", "-rf", final_dir])
    subprocess.call(["mv", "-f", temp_dir, FINAL_ANTISMASH_OUTPUT])

if __name__ == "__main__":
    gbids = []
    gbidfile = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbids.nt.unique.txt"
    ff = open(gbidfile, "r")
    for line in ff:
        gbid = line.strip()
        #if gbid.startswith("AZW"):
       	gbids.append(gbid)

    total = len(gbids)
    print "Total %s" % total
    p = Pool(processes=60)
    #p.map(f, gbids)
    for _ in tqdm.tqdm(p.imap_unordered(f, gbids), total=len(gbids)):
        pass




"""
    l = {} # dict of long species names with short species names as keys
    s = [] # list of short species names
    x = 100000

    for gbid in gbids:
    for genome in os.listdir("/home/gvandova/bacterial_genomes_uncompressed/"):
        #l.append(genome)
        x -= 1
        if x <= 0:
            break
        g = genome.split("_")
        species = g[0] + "_" + g[1]
        if species not in s:
            s.append(species)
            l[species] = genome


    s.sort()
    print ("\n ".join(s))
    print len(l)
    #print len(l.values())
    #p.map(f, l.values())
"""
