#!/usr/bin/env python

# Extract genbank ids with complete genomes

from Bio import SeqIO
import tqdm
import glob
import sys

def get_complete_genomes(outfile, gbfilenames):
    t = open(outfile, "w")
    count = 0

    for gbfilename in tqdm.tqdm(gbfilenames):
        for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
            complete_genome = False
            count += 1
            gbid = record.id.split(".")[0]
            descr = record.description
            if "complete genome" in descr.lower():
                complete_genome = True

            print count, gbid, complete_genome, len(record.seq), descr
            seq = record.seq
            t.write("%s\t%s\t%s" % (gbid, complete_genome, len(record.seq)))
            t.write("\n")
    t.close()




outfile_assemblies = "genome_lengths_assembly.txt"
gbfilenames_assemblies = glob.glob("assembly_gb/*.gbff")

outfile_nt = "genome_lengths_nt.txt"
gbfilenames_nt = glob.glob("gbdir/*.gb")

get_complete_genomes(outfile_assemblies, gbfilenames_assemblies)
get_complete_genomes(outfile_nt, gbfilenames_nt)
