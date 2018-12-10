#!/usr/bin/env python
from Bio import SeqIO
import sys

# To execute script: ./target_ks_tandem.py <filtered_filename> <dist cutoff>
# For example: ./target_ks_tandem.py out.targets.9.filtered 10000

# Max distance of target to KS
DIST_CUTOFF = int(sys.argv[2])

filtered_filename = sys.argv[1]
outfilename = filtered_filename + "." + str(DIST_CUTOFF)
ff = open(outfilename, "w")

targets_to_coord = []
f = open(filtered_filename).readlines()

# FabB/F	KJ189772	2	22654-53136	t2fas	AIW55624.1_ptmP3	1	404	397.0	404.0	404	0.98	0.0
# qseqid, gbid, clusternum, coord, clustertype, protname, prot_start, prot_end, nident, qlen, slen, identity, evalue)

for line in f:
    qseqid, gbid, clusternum, coord, clustertype, protname, prot_start, prot_end, nident, qlen, slen, identity, evalue = line.strip().split("\t")
    if prot_start == "<0":
        prot_start = -1
    prot_start = int(prot_start)
    if ">" in prot_end:
        prot_end = prot_end.split(">")[1]
    prot_end = int(prot_end)
    targets_to_coord.append((qseqid, gbid, clusternum, coord, clustertype, protname, prot_start, prot_end, nident, qlen, slen, identity, evalue))

gbids_to_coord = []

# Find KS coordinates from Blast search
# fasta file from Blast search
fasta_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"

# Head of cdhit file:
# >AVFP01000283.1__724_1992_Microbial_mat_metagenome_scaffold_282__whole_genome_shotgun_sequence_0_1_9914_7e-169
# IAIIGMSGIFPDAEDVQTYWNNLCQGR
# >AM746676___5843905_5845200__0_-1_13033779_0.0

for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
        gbidfull = record.id
        try:
            gbid = gbidfull.split(".")[0]
            coord1 = gbidfull.split("__")[1]
            coord = coord1.split("_")
            start = int(coord[0])
            end = int(coord[1])
            gbids_to_coord.append((gbid, start, end))

        except:
        #if "___" in gbidfull:
            gbid = gbidfull.split("___")[0]
            coord1 = gbidfull.split("___")[1]
            coord2 = coord1.split("__")[0]
            coord = coord2.split("_")
            start = int(coord[0])
            end = int(coord[1])
            gbids_to_coord.append((gbid, start, end))

print len(gbids_to_coord)

for i in xrange(len(targets_to_coord)):
    # print targets_to_coord[i]
    # target_id, gbid1, clusternum, coord, clustertype, protname, target_start,
    # target_end, nident, qlen, slen, identity, evalue = targets_to_coord[i]
    seq_id = targets_to_coord[i][1]
    target_start = int(targets_to_coord[i][6])
    target_end = int(targets_to_coord[i][7])
    # print seq_id, target_start, target_end

    for j in xrange(len(gbids_to_coord)):
        gbid, gb_start, gb_end = gbids_to_coord[j]
        # print gbid
        if seq_id == gbid:
            # Check if targets are withing Xkb of a KS identified from
            # the initial blast search
            dist1 = abs(gb_start - target_end)
            dist2 = abs(target_start - gb_end)
            dist = min(dist1, dist2)
            cluster_start, cluster_end = targets_to_coord[i][3].split("-")
            cluster_start = int(cluster_start)
            cluster_end = int(cluster_end)
            cluster_len = abs(cluster_start - cluster_end)
            if dist < DIST_CUTOFF:
                print seq_id, DIST_CUTOFF, dist
                # print seq_id, dist1, gb_start, gb_end, dist
                ff.write(str(cluster_len))
                ff.write("\t")
                ff.write("\t".join(map(str, targets_to_coord[i])))
                ff.write("\t%d\t%d\t%d" % (gb_start, gb_end, dist))
                ff.write("\n")
                # print "<10kb", target_id, abs(gb_start - target_end), \ gb_start, gb_end

ff.close()
