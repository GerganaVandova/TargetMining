#!/usr/bin/env python
from Bio import SeqIO
import js2py
from js2py import EvalJs
import os
import sys
import tqdm
from multiprocessing import Pool
from collections import defaultdict

gbids_to_coord = defaultdict(list)

fasta_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
# Head of cdhit file:
# >AVFP01000283.1__724_1992_Microbial_mat_metagenome_scaffold_282__whole_genome_shotgun_sequence_0_1_9914_7e-169
# IAIIGMSGIFPDAEDVQTYWNNLCQGR
# >AM746676___5843905_5845200__0_-1_13033779_0.0

count = 0

context = EvalJs()

# Output file with number of gbids on which antismash was run
outfilefaa = "sequences.faa.21k.coord.fasta"
ff = open(outfilefaa, "w")

antismash_dir = "antismash_output_assemblies_all/"
for gbidfull in tqdm.tqdm(os.listdir(antismash_dir)):
    #if not gbidfull.startswith("CP012600"):
    #    continue
    #print gbidfull

    if len(gbidfull.split("_")) > 2:
        gbid = gbidfull.rsplit("_", 2)[0]
        blast_coord_start, blast_coord_end = gbidfull.split("_")[1:]
        #print "here", gbid, blast_coord_start, blast_coord_end
    else:
        gbid = gbidfull
    
    #if gbid != "CP012600":
    #    continue
    
    #print gbid
    count += 1
    print "Parsing %s Number %s: " % (gbid, count)
    filename = antismash_dir + gbidfull + "/geneclusters.js"
    #print filename

    if not os.path.exists(filename):
        continue

    f = open(filename, 'r')
    data_js = f.read()
    context.execute(data_js)
    geneclusters = context.geneclusters.to_dict()
    details_data = context.details_data.to_dict()


    domains = ["PKS_KS", "PKS_KR", "PKS_DH", "PKS_ER", "KS_AT", "KS_ACP",
               "AMP-binding", "Condensation_LCL", "Condensation_DCL", "Condensation_Starter", "PCP",]

    for cluster_id in details_data.iterkeys():
        print "cluster-id from details_data:", cluster_id
        geneclusterorfs = details_data[cluster_id]["orfs"]
        domain_counts = defaultdict(int)
        for orf in geneclusterorfs:
            for domain in orf["domains"]:
                domain_type = domain["type"]
                for d in domains:
                    if d in domain_type:
                        domain_counts[d] += 1

        cluster_key = "%s|%s" % (gbidfull, cluster_id)
        print gbidfull, cluster_id, cluster_key, domain_counts

        cluster_to_domains[cluster_key] = domain_counts

    
    for cluster_id in geneclusters.iterkeys():
        cluster_key = "%s|%s" % (gbidfull, cluster_id)
        print "%s %s" % (cluster_key, json.dumps(cluster_to_domains.get(cluster_key)))

        # write protein sequences in a fasta file
        ff.write("%s\t%s" % (cluster_key, json.dumps(cluster_to_domains.get(cluster_key))))
        ff.write("\n")
        ff.flush()




    for cluster_id in geneclusters.iterkeys():
        
        print "cluster-id from details_data:", cluster_id
        geneclusterorfs = details_data[cluster_id]["orfs"]
        domain_counts = defaultdict(int)
        for orf in geneclusterorfs:
            for domain in orf["domains"]:
                domain_type = domain["type"]
                for d in domains:
                    if d in domain_type:
                        domain_counts[d] += 1
            
        cluster_key = "%s|%s" % (gbidfull, cluster_id)
        print gbidfull, cluster_id, cluster_key, domain_counts
            
        cluster_to_domains[cluster_key] = domain_counts

        geneclustertype = geneclusters[cluster_id]["type"]
        a_start = int(geneclusters[cluster_id]["start"])
        a_end = int(geneclusters[cluster_id]["end"])
        print "CLuster ", cluster_id, "....."

        if len(gbidfull.split("_")) > 2:
            start = int(blast_coord_start) + a_start
            end = int(blast_coord_start) + a_end
        else:
            start = a_start
            end = a_end
        
        for orfs in geneclusters[cluster_id]["orfs"]:
            a_prot_start = int(orfs["start"])
            a_prot_end = int(orfs["end"])

            if len(gbidfull.split("_")) > 2:
                prot_start = int(blast_coord_start) + a_prot_start
                prot_end = int(blast_coord_start) + a_prot_end
            else:
                prot_start = a_prot_start
                prot_end = a_prot_end

            locus_tag = orfs["locus_tag"]
            description = orfs["description"]
            name1 = description.split("</span><br>")[0]
            name = name1.split("<span class=\"svgene-tooltip-bold\">")[1]
            d = description.split("QUERY=")[1]
            sequence1 = d.split("_LOC=protein")[0]
            sequence = sequence1.split("&LINK")[0]
            # write protein sequences in a fasta file
            ff.write(">%s_%s_%s-%s_%s_%s-%s_%s_%s" % (gbid, cluster_id, 
                                                start,
                                                end,
                                                geneclusters[cluster_id]["type"],
                                                prot_start,                                                                                                                                                                  prot_end,
                                                locus_tag, name))
            ff.write("\n")
            ff.write(sequence)
            ff.write("\n")







