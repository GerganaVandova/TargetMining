#!/usr/bin/env python

# Extract taxa info for each genbank id

from Bio import SeqIO
import glob

# taxafile = "taxa_nt.txt"
# taxafile = "taxa_assembly.txt"
# t = open(taxafile, "w")

#taxafile = "species_nt.txt"
taxafile = "species_assembly.txt"
t = open(taxafile, "w")

count = 0

#gbfilenames = glob.glob("gbdir/*.gb")
gbfilenames = glob.glob("assembly_gb/*.gbff")
# gbfilenames = glob.glob("assembly_gb/LQPQ01000178.gbff")
for gbfilename in gbfilenames:
    # print gbfilename
    f = open(gbfilename, 'r')
    data = f.read()
    try:
        for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
            count += 1
            gbid = record.id.split(".")[0]
            seq = record.seq
            name = str(record.annotations["source"])
            print "%s\t%s\t%s" % (count, gbid, name)
            # taxa = list(record.annotations["taxonomy"])
            # for t in taxa:
            #     print "%s\t" % t

            t.write("%s\t%s" % (gbid, name))
            t.write("\n")

    except:
        print "no record for %s" % gbid
        continue
#t.close()
