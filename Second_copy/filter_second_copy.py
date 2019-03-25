#!/usr/bin/python
import sys
import json
from Bio import SeqIO
import os
import subprocess
from collections import defaultdict

# To run script, specify:
# blastdir - this is where blast results are stored
# genome_lengths_file - tab-deliminated file with genbank id,
#                       T/F for complete genome, and genome size
# antismash_filename - output file of ks_target tandem search
# outf - output file with genbank id, target, number of copies and other params
#  ./filter_second_copy.py > second_copy.filtered.genome

def main():

    PIDENTCUTOFF = 0.3
    PIDENTCUTOFF_FAB = 0.6

    blastdir = "blast_out/"

    pairs = defaultdict(int)
    pairs_to_feats = defaultdict(list)

    # AdmT_ACC.1e-08.KB894406.out
    # AdmT_ACC  Subject_1   48.02   205138  205893  121 304 295743  4e-70
    for blastoutfile in os.listdir(blastdir):
        f = os.path.join(blastdir, blastoutfile)
        targetid, _, gbid = f.split(".")[:3]
        lines = open(f).readlines()
        for line in lines:
            targetid, _, sstart, send, \
                nident, qlen, slen, evalue = line.strip().split("\t")
            pidentcutoff = PIDENTCUTOFF
            if "FAB" in targetid:  # for all fabs from the 92-targets list
            # if targetid == "PtmP3_FabB-F":  # for the FabB-F from the 12-targets list
                pidentcutoff = PIDENTCUTOFF_FAB
            pident = float(nident)/int(qlen)
            if float(pident) < pidentcutoff:
                continue
            pairs[(targetid, gbid)] += 1
            pairs_to_feats[(targetid, gbid)] = [sstart, send,
                                                nident, qlen, slen, evalue]

    # Read genome_lengths file
    # A70054    False   1223
    gbid_to_len = {}
    genome_lengths_file = "../Genbank/genome_lengths.txt"
    lines = open(genome_lengths_file, "r").readlines()
    for line in lines:
        gbid, complete_genome, length = line.strip().split("\t")
        gbid_to_len[gbid] = ((complete_genome, length))

    outf = open("out.second_copy.12.20kb.filtered", "w")

    # Read antismash output file
    # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
    antismash_filename = "../Antismash_gbids/out.12.filtered.20kb"
    antismash_file = open(antismash_filename).readlines()
    gbid_to_target = defaultdict(list)
    for line in antismash_file:
        gbid, targetid = line.strip().split("|")[:2]
        complete_genome, length = gbid_to_len[gbid]
        num_copies = pairs[(targetid, gbid)]
        feats = "\t".join(pairs_to_feats[(targetid, gbid)])
        print "%s\t%s\t%s\t%s\t%s\t%s" % \
            (targetid,
             gbid,
             num_copies,
             length,
             complete_genome,
             feats)

        outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                   (targetid,
                    gbid,
                    num_copies,
                    length,
                    complete_genome,
                    feats))

    outf.close()


if __name__ == "__main__":
    main()
