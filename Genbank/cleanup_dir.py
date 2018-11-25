#!/usr/bin/python2.7
import os
import subprocess
import shutil

gbids_file = "gbids.nt.unique.txt"
gbdir = "gbdir"


gbids = set()
with open(gbids_file, "r") as f:
    for gbid in f.readlines():
        gbids.add(gbid.strip())

print "loaded %s gbids" % len(gbids)

gb_files = os.listdir(gbdir)

print "found %s .gb fieles in %s" % (len(gb_files), gbdir)

for gb_file in gb_files:
    if not gb_file.endswith(".gb"):
        print "Not a .gb file found %s" % gb_file
        continue

    gbid = gb_file.split(".gb")[0]
    if gbid not in gbids:
        print "DELETING %s" % gbid
        os.remove(os.path.join(gbdir, gb_file))


antismash_dir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output"

antismash_files = os.listdir(antismash_dir)

for antismash_file in antismash_files:
    if antismash_file not in gbids:
        print "DELETING %s antismash folders" % antismash_file
        shutil.rmtree(os.path.join(antismash_dir, antismash_file))
