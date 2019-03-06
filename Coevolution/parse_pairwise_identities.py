#!/usr/bin/python
from collections import defaultdict
import sys

def get_targets(antismash_outfilename):
    # Get list of target gene names

    # Head of file:
    # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
    targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target = line.split("|")[1]
        targets.add(target)
    return(targets)


def parse_gbids(gbidfull):
    parts = gbidfull.split("|")
    gbid, target = parts[:2]
    id_short = "|".join([gbid, target])
    return gbid, target, id_short


def parse_blast(target, blast_file):
    # Parse blast results; return pairwise target/ks ids and identity
    pairs = defaultdict(float)
    lines = open(blast_file).readlines()
    for line in lines:
        # print line
        # CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
        # 1	421	421	421	421	0.0
        qseqid, sseqid, sstart, send, nident, qlen, slen, evalue = \
            line.strip().split("\t")
        # print qseqid, sseqid, sstart, send, nident, qlen, slen, evalue

        qgb_id, target_query, qseqid_short = parse_gbids(qseqid)
        sgb_id, target_subject, sseqid_short = parse_gbids(sseqid)

        if target_query != target or target_subject != target:
            # print target, target_query, target_subject
            continue
        nident = float(nident)
        qlen = float(qlen)
        identity = float(nident)/qlen*100

        if qseqid > sseqid:
            continue
        if (qseqid, sseqid) in pairs.keys():
            if identity > pairs[(qseqid, sseqid)]:
                pairs[(qseqid, sseqid)] = identity
            else:
                continue

        pairs[(qseqid, sseqid)] = identity
        # pairs[(qseqid, sseqid)] = (identity, line) # if I want to print blast params

    return pairs


def main():

    # Need a file to read target names, use the 10kb cutoff
    antismash_outfilename = "../Antismash_gbids/out.92.filtered.5kb"

    targets = get_targets(antismash_outfilename)
    # targets = ["borI_Thr-tRNA-syn", "mupM_Ile-tRNA-syn", "PtmP3_FabB-F",
    # "rubR1_TIF", "AdmT_ACC", "SalI_beta_proteasome", "BatG_FabI"]

    print targets

    blast_file_ks = "KS.92.5kb.out"
    blast_file_target = "targets.92.5kb.out"
    outfile = "pairwise_identities.92.5kb.out.short"

    # Mibig set
    # blast_file_ks = "KS.mibig.out"
    # blast_file_target = "targets.mibig.out"
    # outfile = "pairwise_identities.mibig.out"

    f = open(outfile, "w")

    for target in targets:
        pairs_ks = parse_blast(target, blast_file_ks)
        pairs_target = parse_blast(target, blast_file_target)

        for (qseqid_ks, sseqid_ks) in sorted(pairs_ks.keys()):
            qgbid, target_q, qseqid_ks_short = parse_gbids(qseqid_ks)
            sgbid, target_s, sseqid_ks_short = parse_gbids(sseqid_ks)

            qks_start, qks_end = qseqid_ks.split("|")[2:4]
            sks_start, sks_end = sseqid_ks.split("|")[2:4]
            len_qks = int(qks_end) - int(qks_start)
            len_sks = int(sks_end) - int(sks_start)

            distance = int(qseqid_ks.split("|")[-1])

            x = pairs_ks[(qseqid_ks, sseqid_ks)]
            y = pairs_target[(qseqid_ks, sseqid_ks)]
            d = abs(x-y)

            # don't write pairwise idientites of short, misannotated KSs
            if len_qks < 300 or len_sks < 300:
                continue

            # print only points on diagonal:
            # if d > 20:
            #     continue
            #
            # if distance < 5000:
            #     continue
            # To print with tab between ids to check low seq identity scores
            # f.write("%s\t%s\t%s\t%s\t%f\t%f\t%s\n" %
            #         (qseqid_ks,
            #          sseqid_ks,
            #          len_qks,
            #          len_sks,
            #          pairs_ks[(qseqid_ks, sseqid_ks)],
            #          pairs_target[(qseqid_ks, sseqid_ks)],
            #          d))
            #
            # print "%s-%s\t%s\t%s\t%.2f\t%.2f\t%.2f" % \
            #     (qgbid,
            #      sgbid,
            #      len_qks,
            #      len_sks,
            #      pairs_ks[(qseqid_ks, sseqid_ks)],
            #      pairs_target[(qseqid_ks, sseqid_ks)],
            #      d)

            # To print short ids:
            print "%s\t%s\t%s\t%s\t%.2f\t%.2f" % \
                (qseqid_ks_short,
                 sseqid_ks_short,
                 len_qks,
                 len_sks,
                 pairs_ks[(qseqid_ks, sseqid_ks)],
                 pairs_target[(qseqid_ks, sseqid_ks)])

            f.write("%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\n" %
                    (qseqid_ks_short,
                     sseqid_ks_short,
                     len_qks,
                     len_sks,
                     pairs_ks[(qseqid_ks, sseqid_ks)],
                     pairs_target[(qseqid_ks, sseqid_ks)],
                     d))

    f.close()


if __name__ == "__main__":
    main()
