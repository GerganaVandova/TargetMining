#!/usr/bin/env python

# Output a file with a genbank id and the corresponding  assembly gb filename

from Bio import SeqIO
import os

gbids_file = "gbids.assembly.unique.txt"
gbdir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Assemblies/assemblies_gbdir/"


assembly_gbids = set()
with open(gbids_file, "r") as f:
    for assembly_gbid in f.readlines():
        assembly_gbids.add(assembly_gbid.strip())

gbid_to_filename = {}

gb_files = os.listdir(gbdir)
print "found %s .gb fieles in %s" % (len(gb_files), gbdir)

outf = open("assemblyids_filenames.2f", 'a')

for gb_file in gb_files:
    #print gb_file
    if not gb_file.endswith(".gbff"):
        print "Not a .gb file found %s" % gb_file
        continue
    
    for record in SeqIO.parse(open(os.path.join(gbdir, gb_file), "rU"), "genbank"):
        gbid = record.id.split(".")[0]
        gbid_to_filename[gbid] = gb_file

        if gbid in assembly_gbids:
            print "%s\t%s" % (gbid, gbid_to_filename[gbid])
            outf.write(gbid)
            outf.write("\t")
            outf.write(gbid_to_filename[gbid])
            outf.write("\n")

outf.close()

    
