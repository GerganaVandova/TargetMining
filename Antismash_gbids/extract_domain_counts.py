#!/usr/bin/env python
from Bio import SeqIO
import js2py
import json
from js2py import EvalJs
import os
import sys
import tqdm
from multiprocessing import Pool
from collections import defaultdict

count = 0
context = EvalJs()

# Output file with number of gbids on which antismash was run
outfiledomains = "domain_counts_assemblies.txt"
ff = open(outfiledomains, "w")


antismash_dir = "antismash_output_assemblies_all/"
for gbidfull in tqdm.tqdm(os.listdir(antismash_dir)):

    cluster_to_domains = {}

    #if not gbidfull.startswith("CP012600"):
    #    continue
    #print gbidfull

    if len(gbidfull.split("_")) > 2:
        gbid = gbidfull.rsplit("_", 2)[0]
        blast_coord_start, blast_coord_end = gbidfull.split("_")[1:]
        #print "here", gbid, blast_coord_start, blast_coord_end
    else:
        gbid = gbidfull
    
    count += 1
    print "\n\nParsing %s Number %s: " % (gbidfull, count)
    filename = antismash_dir + gbidfull + "/geneclusters.js"
    #print filename

    if not os.path.exists(filename):
        continue

    f = open(filename, 'r')
    data_js = f.read()
    context.execute(data_js)
    geneclusters = context.geneclusters.to_dict()
    details_data = context.details_data.to_dict()

    
    print "len details_Data: ", len(details_data.keys())

    domains = [
        "PKS_KS",
        "PKS_KR",
        "PKS_DH",
        "PKS_ER",
        "KS_AT",
        "KS_ACP",
        "AMP-binding",
        "Condensation_LCL",
        "Condensation_DCL",
        "Condensation_Starter",
        "PCP",
    ]


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

ff.close()

