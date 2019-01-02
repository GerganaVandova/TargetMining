#!/usr/bin/env python

# Extract genbank ids with complete genomes

from Bio import SeqIO
import glob


def get_complete_genomes(outfile, gbfilenames):
    
    t = open(outfile, "w")
    
    complete_genomes = []
    count = 0

    for gbfilename in gbfilenames:
        f = open(gbfilename, 'r')
        data = f.readlines()
        for line in data:
                if "DEFINITION" in line:
                    if "Complete genome" in line or "complete genome" in line or "complete chromosome" in line or "Complete chromosome" in line:
                        print gbfilename, line
                        complete_genomes.append(gbfilename)

    for gbfilename in complete_genomes:
        for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
            count += 1
            gbid = record.id.split(".")[0]
            print count, gbid, len(record.seq), record.annotations["keywords"]
            seq = record.seq
            t.write("%s\t%s" % (gbid, len(record.seq)))
            t.write("\n")
    t.close()


outfile_assemblies  = "complete_genomes_assembly.txt"
gbfilenames_assemblies = glob.glob("assembly_gb/*.gbff")

outfile_nt  = "complete_genomes_nt.txt"
gbfilenames_nt = glob.glob("gbdir/*.gb")

get_complete_genomes(outfile_assemblies, gbfilenames_assemblies)
get_complete_genomes(outfile_nt, gbfilenames_nt)

