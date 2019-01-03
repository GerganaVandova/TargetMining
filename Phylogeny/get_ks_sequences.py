#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict
import tqdm


# Get coordinates of KS from Blast results 
fasta_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
fastagbids_to_coord = defaultdict(list)

# CENS01067892.1__80_1441_marine
# CP000510___1361206_1362423__0_1_4559598_

for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
    gbidfull = record.id
    #print gbidfull
    seq = record.seq
    gbid1, rest = gbidfull.split("__", 1)
    gbid = gbid1.split(".")[0]
    parts = filter(lambda x: x, rest.split('_'))
    start = int(parts[0])
    end = int(parts[1])
    print gbid, start, end
    fastagbids_to_coord[gbid].append((start, end, seq))


# Get coordinates of KS from Antismash_output results 
outfilename = "../Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"
antismashgbids_to_coord = defaultdict(list)

# ['Target', 'Cluster', 'Clusternum', 'Clustercoord', 'Type', 'Gene', 'sstart', 'send', 'nident', 'querylen', 'slen', 'pident', 'evalue', 'KS start', 'KS end', 'Distance', 'Cluster len', 'PKS_KS', 'PKS_KR', 'PKS_DH', 'PKS_ER', 'KS_AT', 'KS_ACP', 'AMP-binding', 'Condensation_LCL', 'Condensation_DCL', 'Condensation_Starter', 'PCP']
# AdmT_ACC  ACXX02000001    1   14512-116691    transatpks-nrps Cpap_3701_acetyl-CoA    31720   32586   118.0   304.0   288 0.39    1e-67   38012   39229   5426102179  13  9   7   0

outfile = open(outfilename).readlines()[1:]
for line in outfile:
    line = line.strip()
    cluster = line.split("\t")[1]
    #if cluster == "OFLW01000124":
    #    print cluster, "found"
    
    ks_start = int(line.split("\t")[13])
    ks_end = int(line.split("\t")[14])
    ks_count = int(line.split("\t")[17])
    #print cluster, ks_start, ks_end, ks_count
    antismashgbids_to_coord[cluster].append((ks_start, ks_end, ks_count))


f = open("out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.fasta", "w")

for gbid in tqdm.tqdm(sorted(antismashgbids_to_coord.keys())):
    #if gbid != "OFLW01000124":
    #    continue
    aks_start, aks_end, ks_count = antismashgbids_to_coord[gbid][0] 
#    print gbid, aks_start, aks_end, ks_count

    ksdomains = fastagbids_to_coord[gbid]
    for ks in ksdomains:
        print ks
        ks_start, ks_end, seq = ks
        print ks_start, ks_end, seq[:10]

        if aks_start == ks_start and aks_end == ks_end:
            print ">%s_%s-%s" % (gbid, ks_start, ks_end)
            print seq[:10]
            f.write(">%s_%s-%s" % (gbid, ks_start, ks_end))
            f.write("\n")
            f.write(str(seq))
            f.write("\n")
f.close()


