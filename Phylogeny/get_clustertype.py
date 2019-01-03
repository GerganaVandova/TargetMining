#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
import tqdm

gbid_to_feat = defaultdict(list)

# Get cluster type and target  from antismash output file
# AdmT_ACC  ACXX02000001    1   14512-116691    transatpks-nrps Cpap_3701_acetyl-CoA    31720   32586   118.0   304.0   288 0.39    1e-67   38012   39229   5426102179

antismashfilename = "../Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"
antismashfile = open(antismashfilename).readlines()[1:]

for line in antismashfile:
    feats = line.split("\t")
    target, gbid, cluster_num, coord, cluster_type = feats[:5]
    ks_start = feats[13]
    ks_end = feats[14]
    gbidfull = gbid + "_" + ks_start + "-" + ks_end
    gbid_to_feat[gbidfull] = ((target, cluster_type))

#f = open("out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.cluster_type", "w")
f = open("out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.descr", "w")

# Read KS fasta file and write cluster type and target for each gbid
# >ACXX02000001_38012-39229

fasta_file = "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.fasta"

for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
    gbidfull = record.id
    feats = gbidfull+ "_" + "_".join(gbid_to_feat[gbidfull])
    print feats
    #f.write("%s\t%s\n" % (gbidfull, "_".join(gbid_to_feat[gbidfull])))
    f.write("%s\t%s\n" % (gbidfull, feats))

f.close()

