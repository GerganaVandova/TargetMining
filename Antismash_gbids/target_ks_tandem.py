#!/usr/bin/env python
from Bio import SeqIO
import sys
from collections import defaultdict


def parse_fasta(fasta_file):

    name_to_seq = {}
    with open(fasta_file, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        for record in records:
            # >QEUW01000234|cluster-1|1-12447|t1pks-nrps|6149-6570|DCC55_25085|PKS_KS
            name = record.id
            gbid, cluster_num, cluster_coord, cluster_type, \
                gene_coords, gene_name = name.split("|", 5)
            gene_start, gene_end = gene_coords.split("-")

            seq = record.seq
            name_to_seq[(gbid, gene_start, gene_end)] = seq
            print name, seq[:10]
    return name_to_seq


def main():
    DIST_CUTOFF = 20000
    cluster_genes_file = "cluster_genes.all.fasta"
    ks_file = "ks.all.fasta"
    target_blast_file = "out.12.filtered"

    f1 = open("KS.12.20kb.fasta", "w")
    f2 = open("targets.12.20kb.fasta", "w")
    f3 = open("out.12.filtered.20kb", "w")
    f4 = open("out.12.filtered.20kb.noks", "w")

    min_distance = {}
    data = {}

    with open(ks_file, "rU") as handle:
        ks_records = list(SeqIO.parse(handle, "fasta"))

    with open(target_blast_file, "rU") as handle:
        target_blasts = handle.readlines()

    targets = []
    for line in target_blasts:
        target_name, info, pident, evalue = line.split("\t")
        gbid, cluster_num, cluster_coord, cluster_type, \
            target_coords, gene_name = info.split("|")

        # if gbid != "FRDC01002724": # added for debugging purposes
        #     continue
        target_start, target_end = map(int, target_coords.split("-"))
        targets.append((gbid, target_start, target_end, line))

    print targets # added for debugging purposes

    for target in targets:
        t_gbid, target_start, target_end, target_info = target
        target_n = target_info.split("\t")[0]

        found = False
        for ks_record in ks_records:
            name = ks_record.id
            gbid, cluster_num, cluster_coord, cluster_type, \
                ks_coords, gene_name, domain_type = name.split("|")
            ks_start, ks_end = map(int, ks_coords.split("-"))

            if gbid != t_gbid:
                continue

            dist1 = abs(ks_start - target_end)
            dist2 = abs(target_start - ks_end)
            dist = min(dist1, dist2)

            found = True

            if dist > DIST_CUTOFF:
                continue

            key = (gbid, cluster_coord)
            if min_distance.get(key) is None or min_distance.get(key) > dist:
                min_distance[key] = dist
                # data[key] = "|".join(map(str, [gbid, ks_start, ks_end, target_start, target_end, dist]))
                data[key] = "|".join(map(str, [gbid,
                                               target_n,
                                               ks_start,
                                               ks_end,
                                               target_start,
                                               target_end,
                                               cluster_num,
                                               cluster_type,
                                               cluster_coord,
                                               dist
                                               ]))

        if not found:
            print "no KS found for target %s" % str(target)
            f4.write("%s\n" % str(target))

    ks_to_seq = parse_fasta(ks_file)
    target_to_seq = parse_fasta(cluster_genes_file)

    for v in data.itervalues():
        f3.write("%s\n" % v)
        print v  # > in target_ks_tandem.morethan50kb.out
        # CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
        gbid, target_name, ks_start, ks_end, target_start, target_end = v.split("|")[:6]
        # print gbid, target_name, ks_start, ks_end, target_start, target_end
        print ">%s\n%s" % (v, ks_to_seq[(gbid, ks_start, ks_end)])
        print ">%s\n%s" % (v, target_to_seq[(gbid, target_start, target_end)])
        f1.write(">%s\n%s\n" % (v, ks_to_seq[(gbid, ks_start, ks_end)]))
        f2.write(">%s\n%s\n" % (v, target_to_seq[(gbid, target_start, target_end)]))

    f1.close()
    f2.close()
    f3.close()
    f4.close()

if __name__ == "__main__":
    main()
