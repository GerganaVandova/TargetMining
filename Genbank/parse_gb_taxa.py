#!/usr/bin/env python

# Extract taxa info for each genbank id

from Bio import SeqIO
import glob

#taxafile = "taxa_nt.txt"
#taxafile = "taxa_assembly.txt"
#t = open(taxafile, "w")

count = 0

gbfilenames = glob.glob("gbdir/*.gb")
#gbfilenames = glob.glob("assembly_gb/*.gbff")
# gbfilenames = glob.glob("gbdir/NZ_KE354369.gb")
for gbfilename in gbfilenames:
    # print gbfilename
    f = open(gbfilename, 'r')
    data = f.read()
    try:
        for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
            count += 1
            gbid = record.id.split(".")[0]
            print "\n%s\t" % gbid,
            seq = record.seq
            taxa = list(record.annotations["taxonomy"])
            for t in taxa:
                print "%s\t" % t,
            
            #t.write("%s\t%s" % (gbid, str(record.annotations["taxonomy"])))
            #t.write("\n")
            # print gbid
            # print record.annotations["taxonomy"]
    except:
        continue
#t.close()
