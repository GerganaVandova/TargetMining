#!/usr/bin/env python
from Bio import SeqIO
import js2py
from js2py import EvalJs
import os
import sys
import tqdm
from multiprocessing import Pool



gbids_to_coord = []

fasta_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
# Head of cdhit file:
# >AVFP01000283.1__724_1992_Microbial_mat_metagenome_scaffold_282__whole_genome_shotgun_sequence_0_1_9914_7e-169
# IAIIGMSGIFPDAEDVQTYWNNLCQGR
# >AM746676___5843905_5845200__0_-1_13033779_0.0


for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
        gbidfull = record.id
        try:
            gbid = gbidfull.split(".")[0]
            coord1 = gbidfull.split("__")[1]
            coord = coord1.split("_")
            start = int(coord[0])
            end = int(coord[1])
            gbids_to_coord.append((gbid, start, end))

        except:
        #if "___" in gbidfull:
            gbid = gbidfull.split("___")[0]
            coord1 = gbidfull.split("___")[1]
            coord2 = coord1.split("__")[0]
            coord = coord2.split("_")
            start = int(coord[0])
            end = int(coord[1])
            gbids_to_coord.append((gbid, start, end))
        
#for i in gbids_to_coord:
#    print i

count = 0

context = EvalJs()

# Output file with number of gbids on which antismash was run
outfilefaa = "sequences.faa.89k.coord.fasta"
ff = open(outfilefaa, "w")

antismash_dir = "antismash_output/"
for gbid in os.listdir(antismash_dir):
    count += 1
    print "Parsing %s Number %s: " % (gbid, count)
    filename = antismash_dir + gbid + "/geneclusters.js"
    # print filename

    if not os.path.exists(filename):
        continue

    f = open(filename, 'r')
    data_js = f.read()
    context.execute(data_js)
    geneclusters = context.geneclusters.to_dict()
    details_data = context.details_data.to_dict()

    for cluster_id in geneclusters.iterkeys():
        geneclustertype = geneclusters[cluster_id]["type"]
        start = int(geneclusters[cluster_id]["start"])
        end = int(geneclusters[cluster_id]["end"])
        # print "CLuster ", cluster_id, "....."

        # for j in xrange(len(gbids_to_coord)):
            # print start, abs(start - int(gbids_to_coord[j][1])
            # if gbids_to_coord[j][0] == gbid  and abs(start - int(gbids_to_coord[j][1])) < 100000:
        print gbid, cluster_id, geneclustertype
        for orfs in geneclusters[cluster_id]["orfs"]:
            prot_start = int(orfs["start"])
            prot_end = int(orfs["end"])
            locus_tag = orfs["locus_tag"]
            description = orfs["description"]
            name = description.split("<br>")[0]
            d = description.split("QUERY=")[1]
            sequence1 = d.split("_LOC=protein")[0]
            sequence = sequence1.split("&LINK")[0]
            # print ">%s_%s_%s-%s_%s_%s_%s" % (gbid, cluster_id,
                                            #  geneclusters[cluster_id]["start"],
                                            #  geneclusters[cluster_id]["end"],
                                            #  geneclusters[cluster_id]["type"],
                                            #  locus_tag, name)
            # print sequence

            # write protein sequences in a fasta file
            ff.write(">%s_%s_%s-%s_%s_%s-%s_%s_%s" % (gbid, cluster_id,
                                                geneclusters[cluster_id]["start"],
                                                geneclusters[cluster_id]["end"],
                                                geneclusters[cluster_id]["type"],
                                                prot_start,
                                                prot_end,
                                                locus_tag, name))
            ff.write("\n")
            ff.write(sequence)
            ff.write("\n")







