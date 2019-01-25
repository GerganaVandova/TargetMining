#!/usr/bin/python
import sys
import json

# Append taxa to out file. To execute script:
# ./add_taxa.py ../data/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.unique.276 ../data/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.unique.276.taxa

id_to_features = {}
# /mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/taxa.txt
# NZ_JYJF01001045	['Bacteria', 'Actinobacteria', 'Pseudonocardiales', 'Pseudonocardiaceae', 'Saccharothrix']
# taxa_filename = sys.argv[1]
taxa_filename = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/taxa.txt"
taxa_file = open(taxa_filename).readlines()
for line in taxa_file:
    line = line.strip()
    #print line
    parts = line.split("\t")
    gbid  = parts[0]
    id_to_features[gbid] = parts[1:]
    #print id_to_features[gbid]


#AdmT_ACC   ACXX02000001    1   14512-116691    transatpks-nrps Cpap_3701_acetyl-CoA    31720   32586   118.0   304.0   288 0.39    1e-67   38012   39229   5426    102179
input_filename = sys.argv[1]
input_file = open(input_filename).readlines()

output_file = input_filename + ".taxa"
outf = open(output_file, "w")

for line in input_file:
    line = line.strip()
    features = line.split("\t")
    cluster_name = features[1]
    #target, cluster_name, cluster_num = features[:4] # for e coli targets
    feats = str('\t'.join(features))
    if cluster_name not in id_to_features.keys():
        outf.write("%s" % features)
        outf.write("\n")
        print "No taxa", cluster_name
    if cluster_name in id_to_features.keys():
        # taxa = id_to_features[cluster_name]
        print type(id_to_features[cluster_name])
        taxa = '\t'.join(id_to_features[cluster_name])
        print taxa
        #sys.exit(0)
        outf.write("%s\t%s" % (feats, taxa))
        outf.write("\n")
        print feats, "\t", taxa

outf.close()
