#!/usr/bin/python
from collections import defaultdict
import sys
import tqdm


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
        data[target_name].append((gbid1, gbid2, d))
    return data


def get_score(target, gbid, pairwise_data):
    distance = 0
    count = 0
    for gbid1, gbid2, d in pairwise_data[target]:
        if gbid1 == gbid or gbid2 == gbid:
            count += 1
            distance += d
    return (distance, count)


def main():
    f = open("coevolution_scores.12.20kb", "w")
    gbid_to_score = defaultdict()

    pairwise_filename = "pairwise_identities.12.20kb.out"
    pairwise_data = load_pairwise_file(pairwise_filename)
    antismash_outfilename = "../Antismash_gbids/out.12.filtered.20kb"
    gbids_and_targets = get_gbids_and_targets(antismash_outfilename)
    for gbid, target in tqdm.tqdm(sorted(gbids_and_targets)):
        distance, count = get_score(target, gbid, pairwise_data)
        if count != 0:
            score = float(distance)/count
            f.write("%s\t%s\t%s\t%s\n" % (gbid, target, count, score))
    f.close()

if __name__ == "__main__":
    main()
