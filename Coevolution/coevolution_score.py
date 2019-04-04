#!/usr/bin/python
from collections import defaultdict
import sys
import tqdm
import json

def get_gbids_and_targets(antismash_outfilename):
    gbids_and_targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        gbid, target, rest = line.split("|", 2)
        gbids_and_targets.add((gbid, target))
    return gbids_and_targets


def load_pairwise_file(pairwise_filename):
    # target -> [(gbid1, gbid2, d)]
    #### target -> {gbid -> (count, distance)}
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
    red = 0
    red_list = []
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

            # Get redundant clusters
            if pident_ks >= 90:
                red += 1
                if gbid1 != gbid:
                    red_list.append((gbid1))
                elif gbid2 != gbid:
                    red_list.append((gbid2))
    return (distance, count, homologs, homologs_list, pident_list, red, red_list)


def main():
    f = open("coevolution_scores.609.10kb", "w")
    gbid_to_score = defaultdict()

    pairwise_filename = "pairwise_identities.609.10kb.out"
    pairwise_data = load_pairwise_file(pairwise_filename)
    antismash_outfilename = "../Antismash_gbids/out.609.filtered.10kb"
    gbids_and_targets = get_gbids_and_targets(antismash_outfilename)
    for gbid, target in tqdm.tqdm(sorted(gbids_and_targets)):
        distance, count, homologs, homologs_list, pident_list, \
            red, red_list = get_score(target, gbid, pairwise_data)
        if count != 0:
            score = float(distance)/count
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                    (gbid, target, count, score,
                     homologs, json.dumps(homologs_list), json.dumps(pident_list),
                     red, json.dumps(red_list)))
    f.close()

if __name__ == "__main__":
    main()
