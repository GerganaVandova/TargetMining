#!/usr/bin/env python

# Convert gb file to fasta file, write dna sequence, as not all genbank files
# have predicted protein sequences.

from Bio import SeqIO
import glob

pseqfile = "prot_nt_test.fasta"
fp = open(pseqfile, "w")

dseqfile = "dna_nt_test.fasta"
fd = open(dseqfile, "w")

count_dnaseq = 0
count_protseq = 0
count_noseq = 0
noseq = set()

gbfilenames = glob.glob("gbdir/*.gb")
# gbfilenames = glob.glob("gbdir/NZ_KE354369.gb")
for gbfilename in gbfilenames:
    #print gbfilename
    f = open(gbfilename, 'r')
    data = f.read()
    for record in SeqIO.parse(open(gbfilename, "rU"), "genbank"):
        gbid = record.id.split(".")[0]
	print gbid
    # print record.annotations["organism"]
    # sys.exit(0)
	#print record.features
	sequence = record.seq
	if sequence[0:10] == "NNNNNNNNNN":
            #print gbid, " no DNA sequence"

	    for feature in record.features:
		#print "feature location", feature.location
		coord1 = str(feature.location)
		coord = coord1.split("(")[0]
	    	if feature.type == "CDS":
                    count_protseq += 1
                    #print ("count protseq %s\t gbid %s" % (count_protseq, gbid)) 
		    try:
		        #print "feature qualifiers", feature.qualifiers
	    	        prot_id = "".join(feature.qualifiers["protein_id"])
	    	        prots = "".join(feature.qualifiers["translation"])
			count_protseq += 1
                        #print ("protseq%s\t%s%s" % (count_protseq, gbid, str(prots)[:60]))
		        #fp.write(">%s__%s__%s" % (gbid, coord, prot_id))
		        #fp.write("\n")
		        #fp.write("%s" % str(prots))
		        #fp.write("\n")
		    except:
           		count_noseq += 1
			noseq.add(gbid)
                        print ("count no seq %s\t gbid %s" % (count_noseq, gbid))
	else:
            count_dnaseq += 1
	    #print("count dnaseq %s\t%s" % (count_dnaseq, gbid))
       	    #fd.write(">%s" % gbid)
            #fd.write("\n")
            #fd.write(str(sequence))
            #fd.write("\n")

fp.close()
fd.close()

print "Prot seq", count_protseq
print "Dna seq", count_dnaseq
print "no seq", count_noseq, noseq

print noseq
print len(noseq)
