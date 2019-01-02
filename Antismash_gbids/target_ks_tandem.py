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

min_distance = {}
data = {}

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
#    print gbidfull
    gbid, rest = gbidfull.split("__", 1)
    gbid = gbid.split(".")[0]
    parts = filter(lambda x: x, rest.split('_'))
    start = int(parts[0])
    end = int(parts[1])
#    if gbid.startswith("KT362046"):
#        print gbid
#    print gbid, start, end
#    gbids_to_coord[gbid].append((start, end))
    gbids_to_coord.append((gbid, start, end))


#print len(gbids_to_coord)

def get_key(arr):
    return (arr[0], arr[1], arr[2], arr[3], arr[5])

for i in xrange(len(targets_to_coord)):
    # print targets_to_coord[i]
    # target_id, gbid1, clusternum, coord, clustertype, protname, target_start,
    # target_end, nident, qlen, slen, identity, evalue = targets_to_coord[i]
    seq_id = targets_to_coord[i][1]
    target_start = int(targets_to_coord[i][6])
    target_end = int(targets_to_coord[i][7])
    
    #if seq_id != "KT362046":
    #    continue

    #print seq_id, target_start, target_end

    for j in xrange(len(gbids_to_coord)):
        gbid, gb_start, gb_end = gbids_to_coord[j]
        if seq_id == gbid:
            print targets_to_coord[i]
            print "KS coord: ", gb_start, gb_end
            print "target coords: ", targets_to_coord[i][0], target_start, target_end 
            # Check if targets are within Xkb of a KS identified from
            # the initial blast search
            dist1 = abs(gb_start - target_end)
            dist2 = abs(target_start - gb_end)
            dist = min(dist1, dist2)
            cluster_start, cluster_end = targets_to_coord[i][3].split("-")
            cluster_start = int(cluster_start)
            cluster_end = int(cluster_end)
            cluster_len = abs(cluster_start - cluster_end)
            print "dist", dist, "\n"
            if dist < DIST_CUTOFF:
                print seq_id, DIST_CUTOFF, dist
                # print seq_id, dist1, gb_start, gb_end, dist
                line = "%s\t%s\t%s" % (
                    "\t".join(map(str, targets_to_coord[i])), 
                    "\t".join(map(str, [gb_start, gb_end, dist])),
                    cluster_len)
#                ff.write(line)
#                ff.write("\n")
                key = get_key(targets_to_coord[i])
                if min_distance.get(key) is None or min_distance.get(key) > dist:
                    min_distance[key] = dist
                    data[key] = line 


for v in data.itervalues():
    ff.write(v)
    ff.write("\n")

ff.close()
