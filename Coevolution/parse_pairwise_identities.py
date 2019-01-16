#!/usr/bin/python
from collections import defaultdict


def get_targets(antismash_outfilename):
    # Get list of target gene names

    targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target = line.split("\t")[0]
        targets.add(target)
    return(targets)


def parse_blast(target, blast_file):
    # Parse blast results; return pairwise target/ks ids and identity
    pairs = defaultdict(float)
    lines = open(blast_file).readlines()
    for line in lines:
        # AdmT_ACC.ACXX02000001.14512-116691.KS.38012-39229
        # AdmT_ACC.ACXX02000001.14512-116691.KS.38012-39229
        # 1 406	406	406	406	0.0
        qseqid, sseqid, sstart, send, nident, qlen, slen, evalue = \
            line.strip().split("\t")
        # print qseqid, sseqid, sstart, send, nident, qlen, slen, evalue
        target_query = qseqid.split(".")[0]
        if target_query != target:
            continue
        nident = float(nident)
        qlen = float(qlen)
        identity = float(nident)/qlen
        qseqid_short = qseqid.rsplit(".", 3)[0]
        sseqid_short = sseqid.rsplit(".", 3)[0]
        # if (qseqid, sseqid) in pairs.keys() or (sseqid, qseqid) in pairs.keys():
        #     print "%s\t%s\t%.4f already exists" % (qseqid_short, sseqid_short, identity)
        #     # if identity > pairs[(qseqid, sseqid)]:
        #     #     # print "%s\t%s\t%f  previous identity " % (qseqid_short, sseqid_short, pairs[(qseqid, sseqid)])
        #     #     # print "%s\t%s\t%f found higher seq identity " % (qseqid_short, sseqid_short, identity)
        #     #     pairs[(qseqid, sseqid)] = identity
        #     # if identity > pairs[(sseqid, qseqid)]:
        #     #     # print "%s\t%s\t%f  previous identity " % (qseqid_short, sseqid_short, pairs[(qseqid, sseqid)])
        #     #     # print "%s\t%s\t%f found higher seq identity " % (qseqid_short, sseqid_short, identity)
        #     #     pairs[(sseqid, qseqid)] = identity
        # else:
        #     pairs[(qseqid, sseqid)] = identity
        pairs[(qseqid, sseqid)] = identity
    return pairs


antismash_outfilename = \
    "../Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"

targets = get_targets(antismash_outfilename)
blast_file_ks = "all.KS.out"
blast_file_target = "all.targets.out"

outfile = "pairwise_identities.out"
f = open(outfile, "w")

for target in targets:
    pairs_ks = parse_blast(target, blast_file_ks)
    pairs_target = parse_blast(target, blast_file_target)
    for (qseqid_ks, sseqid_ks) in sorted(pairs_ks.keys()):
        qseqid_k = qseqid_ks.rsplit(".", 2)[0]
        sseqid_k = sseqid_ks.rsplit(".", 2)[0]

        qseqid_k_short = qseqid_ks.rsplit(".", 3)[0]
        sseqid_k_short = sseqid_ks.rsplit(".", 3)[0]

        for (qseqid_target, sseqid_target) in sorted(pairs_target.keys()):
            # AdmT_ACC.ACXX02000001.14512-116691.KS.38012-39229
            qseqid_t = qseqid_target.rsplit(".", 1)[0]
            sseqid_t = sseqid_target.rsplit(".", 1)[0]

            qseqid_t_short = qseqid_target.rsplit(".", 2)[0]
            sseqid_t_short = sseqid_target.rsplit(".", 2)[0]

            if (qseqid_k, sseqid_k) == (qseqid_t, sseqid_t):
                f.write("%s-%s\t%.4f\t%s-%s\t%.4f\n" %
                        (qseqid_k_short,
                         sseqid_k_short,
                         pairs_ks[(qseqid_ks, sseqid_ks)],
                         qseqid_t_short,
                         sseqid_t_short,
                         pairs_target[(qseqid_target, sseqid_target)]))

                print "%s-%s\t%.4f\t%s-%s\t%.4f" % \
                    (qseqid_k_short,
                     sseqid_k_short,
                     pairs_ks[(qseqid_ks, sseqid_ks)],
                     qseqid_t_short,
                     sseqid_t_short,
                     pairs_target[(qseqid_target, sseqid_target)])

f.close()
