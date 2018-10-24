#!/usr/bin/env python
import subprocess
import os
import os.path
from Bio import Entrez
from multiprocessing import Pool


def f(gbid):
    print "\n", "***********", "\n", gbid
    antismashfolder = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output/"
    antismashoutdir = os.path.join(antismashfolder + gbid)
    subprocess.call(["mkdir", "-p", antismashoutdir])
    gbfolder = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbdir"
    local_filename = os.path.join(gbfolder, gbid + ".gb")
    geneclusters_filename = os.path.join(antismashoutdir, "geneclusters.js")
    print geneclusters_filename
    if os.path.exists(geneclusters_filename) == True:
        print geneclusters_filename, " exist"
        return
    antismash = "python /home/gvandova/antismash/run_antismash.py" + " --outputfolder " + \
                antismashoutdir + " " + local_filename
    print antismash
    print "start %s" % gbid
    os.system(antismash)


if __name__ == "__main__":

    gbids = []
    gbidfile = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbids.unique.84233.txt"
    ff = open(gbidfile, "r")
    for line in ff:
        gbid = line.strip()
        gbids.append(gbid)

    p = Pool(processes=60)
    p.map(f, gbids)




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
