#!/usr/bin/python
from collections import defaultdict
import sys
import tqdm
import json
from Bio import SeqIO



def get_target_counts(antismash_ksfasta):
    target_counts = {}
    targets = set()
    for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
        gbidfull = record.id
        target = gbidfull.split("|")[1]
        if target in targets:
            target_counts[target] += 1
        else:
            targets.add(target)
            target_counts[target] = 1
    return (targets, target_counts)


def get_gbids_and_targets(antismash_ksfasta):
    gbids_and_targets = set()
    for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
        gbidfull = record.id
        gbid, target, rest = gbidfull.split("|", 2)
        gbids_and_targets.add((gbid, target))
    return gbids_and_targets


def load_pairwise_file(pairwise_filename):
    data = defaultdict(list)
    lines = open(pairwise_filename).readlines()
    for line in lines:
        qseqid, sseqid, len_qks, len_sks, pident_ks, pident_target, d = \
            line.strip().split("\t")
        d = float(d)
        gbid1, target_name, rest = qseqid.split("|", 2)
        gbid2 = sseqid.split("|")[0]
        data[target_name].append((gbid1, gbid2, d, pident_ks))
    return data


def get_score(target, gbid, pairwise_data):
    distance = 0
    count = 0
    homologs = 0
    homologs_list = []
    pident_list = []
    for gbid1, gbid2, d, pident_ks in pairwise_data[target]:
        if gbid1 == gbid or gbid2 == gbid:
            count += 1
            distance += d
            pident_ks = float(pident_ks)

            # Get homologs
            if pident_ks > 80 and pident_ks < 90:
                homologs += 1
                pident_list.append(pident_ks)
                if gbid1 != gbid:
                    homologs_list.append((gbid1))
                elif gbid2 != gbid:
                    homologs_list.append((gbid2))

    return (distance, count, homologs, homologs_list, pident_list)


def main():
    f = open("coevolution_scores.12.10kb.200bpks", "w")
    gbid_to_score = defaultdict()

    pairwise_filename = "pairwise_identities.12.10kb.200bpks.out"
    pairwise_data = load_pairwise_file(pairwise_filename)
    antismash_ksfasta = "../Antismash_gbids/KS.12.10kb.fasta.cdhit.90"
    targets, target_counts = get_target_counts(antismash_ksfasta)
    for target in target_counts.iterkeys():
        print "%s\t%s" % (target, target_counts[target])

    gbids_and_targets = get_gbids_and_targets(antismash_ksfasta)
    for gbid, target in tqdm.tqdm(sorted(gbids_and_targets)):
        distance, count, homologs, homologs_list, pident_list, \
            = get_score(target, gbid, pairwise_data)
        target_ubiquity = target_counts[target]
        if count != 0:
            score = float(distance)/count
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                    (gbid, target, target_ubiquity, score,
                     homologs, json.dumps(homologs_list), json.dumps(pident_list)))
    f.close()

if __name__ == "__main__":
    main()
