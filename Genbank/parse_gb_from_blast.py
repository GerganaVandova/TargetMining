#!/usr/bin/env python

# For each cluster make a fasta file with all proteins 

def split_fasta(infile):
    # Split fasta file into files for each sequence
    # Store seq files in uniprot_seqs/ folder

    with open(infile, "rU") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
        for record in records:
            name = record.name
            #>BGC0000001|c1|1-1083|+|no_locus_tag|protein_methyltransferase|AEK75490.1
            print name
            cluster = name.split("|")[0]
            outfilename = "clusters_fasta/%s.fasta" % cluster
            f = open(outfilename, 'a')
            f.write(">%s\n" % name)
            f.write(str(record.seq))
            f.write("\n")



gbfilenames = glob.glob("gbdir/*.gb")
for gbfilename in gbfilenames:

infile = 'MIBiG_prot_seqs_1.4.fa'
#infile = 'all_uniprot.fa'
split_fasta(infile)


# Convert gb file to fasta file, write dna sequence, as not all genbank files
# have predicted protein sequences.




from Bio import SeqIO
import glob

pseqfile = "prot_nt.fasta"
fp = open(pseqfile, "w")

dseqfile = "dna_nt.fasta"
fd = open(dseqfile, "w")

count_dnaseq = 0
count_protseq = 0
count_noseq = 0
noseq = []

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
	    for feature in record.features:
		#print "feature location", feature.location
		coord1 = str(feature.location)
		coord = coord1.split("(")[0]
	    	if feature.type == "CDS":
		    try:
		        #print "feature qualifiers", feature.qualifiers
	    	        prot_id = "".join(feature.qualifiers["protein_id"])
	    	        prots = "".join(feature.qualifiers["translation"])
			count_protseq += 1
                        print ("protseq%s\t%s%s" % (count_protseq, gbid, str(prots)[:60]))
		        fp.write(">%s__%s__%s" % (gbid, coord, prot_id))
		        fp.write("\n")
		        fp.write("%s" % str(prots))
		        fp.write("\n")
		    except:
           		count_noseq += 1
			noseq.append(gbid)
	else:
            count_dnaseq += 1
	    print("dnaseq%s\t%s\n%s" % (count_dnaseq, gbid, sequence[0:60]))
	    count_dnaseq += 1
       	    fd.write(">%s" % gbid)
            fd.write("\n")
            fd.write(str(sequence))
            fd.write("\n")

fp.close()
fd.close()

print "Prot seq", count_protseq
print "Dna seq", count_dnaseq
print "no seq", count_noseq, noseq
