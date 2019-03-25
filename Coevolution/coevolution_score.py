#!/usr/bin/python
from collections import defaultdict
import sys
import tqdm

#
# def get_targets(antismash_outfilename):
#     # Get list of target gene names
#     # Head of file:
#     # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
#     targets = set()
#     outfile = open(antismash_outfilename).readlines()[1:]
#     for line in outfile:
#         line = line.strip()
#         target = line.split("|")[1]
#         targets.add(target)
#     return(targets)


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
        # print line
        qseqid, sseqid, len_qks, len_sks, pident_ks, pident_target, d = \
            line.strip().split("\t")
        d = float(d)
        gbid1, target_name, rest = qseqid.split("|", 2)
        gbid2 = sseqid.split("|")[0]
        data[target_name].append((gbid1, gbid2, d))
    return data

def get_score(target, gbid, pairwise_data):
    # Parse blast results; return pairwise target/ks ids and identity
    distance = 0
    count = 0
    for gbid1, gbid2, d in pairwise_data[target]:
        if gbid1 == gbid or gbid2 == gbid:
            count += 1
            distance += d
            # print count, distance, gbid1, gbid2, target_name
    return (distance, count)


def main():
    f = open("coevolution_scores.609.10kb", "w")
    # targets = ["PtmP3_FabB-F", "agnB2_Leu-tRNA-syn"]
    # gbids = ["KE354369", "LN879412", "AM889285"]
    # target = "PtmP3_FabB-F"
    # gbid = "KE354369"
    # print get_score(target, gbid, pairwise_filename)

    gbid_to_score = defaultdict()

    pairwise_filename = "pairwise_identities.609.10kb.out"
    pairwise_data = load_pairwise_file(pairwise_filename)
    antismash_outfilename = "../Antismash_gbids/out.609.filtered.10kb"
    gbids_and_targets = get_gbids_and_targets(antismash_outfilename)
    # targets = get_targets(antismash_outfilename)
    for gbid, target in tqdm.tqdm(sorted(gbids_and_targets)):
        # for gbid in sorted(gbids):
        distance, count = get_score(target, gbid, pairwise_data)
        # print gbid, target, distance, count
        if count != 0:
            score = float(distance)/count
            #print gbid, target, count, score
            f.write("%s\t%s\t%s\t%s\n" % (gbid, target, count, score))
#             gbid_to_score[gbid] = target, count, score
    #
    # for gbid in sorted(gbid_to_score.keys()):
    #     print gbid, target, count, score
    f.close()

if __name__ == "__main__":
    main()
