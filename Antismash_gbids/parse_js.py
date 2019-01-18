#!/usr/bin/env python
from Bio import SeqIO
import js2py
from js2py import EvalJs
import os
import sys
import tqdm
from multiprocessing import Pool
from collections import defaultdict


def get_orfs(gbidfull, antismash_dir, gene_outfile, ks_outfile, stdoutfile):
    filename = os.path.join(antismash_dir, gbidfull, "geneclusters.js")

    if not os.path.exists(filename):
        print "Not found %s" % filename
        return

    gbid, gbid_start, gbid_end = parse_gbidfull(gbidfull)

    # Read antismash json file
    f = open(filename, 'r')
    data_js = f.read()
    context = EvalJs()
    context.execute(data_js)
    geneclusters = context.geneclusters.to_dict()
    details_data = context.details_data.to_dict()

    for cluster_id in geneclusters.iterkeys():
        # Parse only pks clusters:
        clustertype = geneclusters[cluster_id]["type"]
        pkss = ['t1pks', 'transatpks']
        if not any(pks in clustertype for pks in pkss):
            continue

        # Correct gene coord if gbid sequence was split
        antismash_start = int(geneclusters[cluster_id]["start"])
        antismash_end = int(geneclusters[cluster_id]["end"])
        cluster_start = gbid_start + antismash_start
        cluster_end = gbid_start + antismash_end

        print "Cluster ", cluster_id, cluster_start, cluster_end, "....."
        locus_coords = {}
        for orfs in geneclusters[cluster_id]["orfs"]:
            antismash_prot_start = int(orfs["start"])
            antismash_prot_end = int(orfs["end"])
            #print antismash_prot_start, antismash_prot_end
            # Correct gene coord if gbid sequence was split
            prot_start = gbid_start + antismash_prot_start
            prot_end = gbid_start + antismash_prot_end

            locus_tag = orfs["locus_tag"]
            locus_coords[locus_tag] = (prot_start, prot_end)
            description = orfs["description"]
            name1 = description.split("</span><br>")[0]
            name = name1.split("<span class=\"svgene-tooltip-bold\">")[1]
            d = description.split("QUERY=")[1]
            sequence1 = d.split("_LOC=protein")[0]
            sequence = sequence1.split("&LINK")[0]

            # Write protein sequences in a fasta file
            fasta_id = ">%s|%s|%s-%s|%s|%s-%s|%s" % \
                (gbid, cluster_id, cluster_start, cluster_end,
                 clustertype, prot_start, prot_end, locus_tag)
            gene_outfile.write("%s\n%s\n" % (fasta_id, sequence))

        if not details_data:
            print "%s\t%s\tEmpty details_data" % \
                (gbidfull, abs(cluster_end - cluster_start))
            stdoutfile.write("%s\t%s\tEmpty details_data\n" %
                             (gbidfull, abs(cluster_end - cluster_start)))
            continue
        for orfs in details_data[cluster_id]["orfs"]:
            locus_tag = orfs["id"]
            for domain in orfs["domains"]:
                domain_type = domain["type"]
                if domain_type != "PKS_KS":
                    continue
                ks_seq = domain["sequence"]

                ks_start = locus_coords[locus_tag][0] + int(domain["start"])
                ks_end = locus_coords[locus_tag][0] + int(domain["end"])
                print "KS %s %s %s" % (ks_start, ks_end, ks_seq[:20])

                # Write protein sequences in a fasta file
                fasta_id = ">%s|%s|%s-%s|%s|%s-%s|%s|%s" % \
                    (gbid, cluster_id, cluster_start, cluster_end,
                     clustertype, ks_start, ks_end, locus_tag, domain_type)
                ks_outfile.write("%s\n%s\n" % (fasta_id, ks_seq))

#
# def read_cluster_genes(outfile_cluster_genes):
#     fasta_to_seq = {}
#     with open(outfile_cluster_genes, "rU") as handle:
#         records = list(SeqIO.parse(handle, "fasta"))
#         for record in records:
#             name = record.id
#             seq = record.seq
#             fasta_to_seq[name] = seq
#     return fasta_to_seq

#
#
# def get_ks_seqs(filename, fasta_to_seq, outfile_ks, gbid):
#
#     # Read antismash json file
#     f = open(filename, 'r')
#     data_js = f.read()
#     context.execute(data_js)
#     geneclusters = context.geneclusters.to_dict()
#     details_data = context.details_data.to_dict()
#
#     for cluster_id in geneclusters.iterkeys():
#         # Parse only pks clusters:
#         clustertype = geneclusters[cluster_id]["type"]
#         pkss = ['t1pks', 'transatpks']
#         if not any(pks in clustertype for pks in pkss):
#             continue
#
#         # Correct gene coord if gbid sequence was split
#         antismash_start = int(geneclusters[cluster_id]["start"])
#         antismash_end = int(geneclusters[cluster_id]["end"])
#         cluster_start = int(abs_coord_start) + antismash_start
#         cluster_end = int(abs_coord_start) + antismash_end
#
#         partial_key = "%s|%s|%s-%s|%s" % (gbid, cluster_id, cluster_start,
#                                    cluster_end, clustertype)
#         for key in fasta_to_seq.keys():
#             # print "key: ", key
#             # print "partial key: ", partial_key
#             if partial_key in key:
#                 print "FOUND MATCH", partial_key, key
#                 orf_seq = fasta_to_seq[key]
#                 orf_coord = key.split("|")[4]
#                 orf_start, orf_end = orf_coord.split("-")
#
#                 for cluster_id in geneclusters.iterkeys():
#                     # Parse only pks clusters:
#                     clustertype = geneclusters[cluster_id]["type"]
#                     pkss = ['t1pks', 'transatpks']
#                     if not any(pks in clustertype for pks in pkss):
#                         continue
#
#                     # Extract KS sequences
#                     orfs = details_data[cluster_id]["orfs"]
#                     for orf in orfs:
#                         for domain in orf["domains"]:
#                             domain_type = domain["type"]
#                             if domain_type != "PKS_KS":
#                                 continue
#                             ks_seq = domain["sequence"]
#                             antismash_ks_start = int(domain["start"])
#                             antismash_ks_end = int(domain["end"])
#                             if ks_seq in orf_seq:
#                                 # Correct gene coord if gbid sequence was split
#                                 ks_start = int(orf_start) + antismash_ks_start
#                                 ks_end = int(orf_start) + antismash_ks_end
#                                 print "antismash_ks_start, antismash_ks_end, ks_strt, ks_end:", antismash_ks_start, antismash_ks_end, ks_start, ks_end
#
#                                 # Write KS sequences in a fasta file
#                                 fasta_ksid = ">%s|%s|%s-%s|%s|%s-%s|PKS_KS" % \
#                                     (gbid, cluster_id, cluster_start, cluster_end,
#                                      clustertype, ks_start, ks_end)
#                                 ffks.write("%s\n%s\n" % (fasta_ksid, ks_seq))

def parse_gbidfull(gbidfull):
    # Check if gbid was split for the antismash step, then the folder name
    # should be [gbid]_[coordstart]_[coordend]: CP012600_1054705_1355970
    # These are the absolute genbank coordinates
    # Check if gbid sequence was split and if so correct cluster coord
    if len(gbidfull.split("_")) > 2:
        gbid = gbidfull.rsplit("_", 2)[0]
        abs_coord_start, abs_coord_end = gbidfull.split("_")[1:]
        return gbid, int(abs_coord_start), int(abs_coord_end)
    else:
        return  gbidfull, 0, None


def main():
    count = 0
    # antismash_dir = "antismash_output/"
    antismash_dir = "antismash_output_assemblies_all/"

    # outfile_cluster_genes = "cluster_genes.89k.fasta"
    outfile_cluster_genes = "cluster_genes.21k.fasta.test"
    gene_out_f = open(outfile_cluster_genes, "w")

    outfile_ks = "ks.21k.fasta.test"
    ks_out_f = open(outfile_ks, "w")

    stdoutfile = "stdout.txt"
    stdout_f = open(stdoutfile, "w")

    # infile_cluster_genes = "cluster_genes.21k.fasta"

    for gbidfull in tqdm.tqdm(os.listdir(antismash_dir)):
        if not gbidfull.startswith("JOIR01000095"):
            continue
        print gbidfull

        count += 1
        print "Parsing %s Number %s: " % (gbidfull, count)
        get_orfs(gbidfull, antismash_dir, gene_out_f, ks_out_f, stdout_f)

if __name__ == "__main__":
    main()
