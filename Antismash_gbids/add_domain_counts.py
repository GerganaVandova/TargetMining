#!/usr/bin/python
import sys
import json

#To execute script: ./add_domain_counts.py all.cluster.counts.txt out.targets.12.eval.1-8.pident.30.filtered.10000.allpks out.targets.12.eval.1-8.pident.30.filtered.20000.allpks.domaincounts

id_to_features = {}
# domain_counts.txt
#CP012600_1054705_1355970|cluster-1  {"Condensation_DCL": 2, "PCP": 7, "KS_AT": 1, "AMP-binding": 6, "PKS_KR": 2, "PKS_KS": 1, "Condensation_LCL": 3, "Condensation_Starter": 1}
#CP012600_1054705_1355970|cluster-2  null
#CP012600_524517_825731|cluster-1    null
#CP012600_2591087_2892367|cluster-1  {"PKS_DH": 1, "PCP": 2, "AMP-binding": 2, "PKS_KR": 1, "PKS_KS": 1, "Condensation_LCL": 1, "KS_AT": 1}
#CP012600_2591087_2892367|cluster-2  null


domains = [
        "PKS_KS",
        "PKS_KR",
        "PKS_DH",
        "PKS_ER",
        "KS_AT",
        "KS_ACP",
        "AMP-binding",
        "Condensation_LCL",
        "Condensation_DCL",
        "Condensation_Starter",
        "PCP",
    ]

ks_count_filename = sys.argv[1]  # domain_counts.txt
ks_count_file = open(ks_count_filename).readlines()
for line in ks_count_file:
    line = line.strip()
    cluster_key, rest  = line.split("\t")
    x = json.loads(rest) or {}
    id_to_features[cluster_key] = x
    

#AdmT_ACC    AM420293    1   8925-62106  t1pks-otherks   ctg1_25_-   45467   47173   125.0   304.0   568 0.41    1e-79   36533   37813   7654    53181
#AdmT_ACC    AM420293    1   8925-62106  t1pks-otherks   SACE_0026_acetyl    45467   47173   125.0   304.0   568 0.41    1e-79   36533   37813   7654    53181
#GyrB-R  CP012600    1   1178914-1248724 t1pks-nrps  ctg1_140_-  1185219 1187168 305.0   677.0   649 0.45    0.0 1204705 1205970 17537   69810

input_filename = sys.argv[2]  # Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.20000.allpks 
input_file = open(input_filename).readlines()

output_file = sys.argv[3]
outf = open(output_file, "w")
outf.write("Target\tCluster\tClusternum\tClustercoord\tType\tGene\tsstart\tsend\tnident\tquerylen\tslen\tpident\tevalue\tKS start\tKS end\tDistance\tCluster len\tPKS_KS\tPKS_KR\tPKS_DH\tPKS_ER\tKS_AT\tKS_ACP\tAMP-binding\tCondensation_LCL\tCondensation_DCL\tCondensation_Starter\tPCP\n")

for line in input_file:
    line = line.strip()
    features = line.split("\t")
    cluster_name, out_cluster_num, coords = features[1:4]
    out_coord_start, out_coord_end = coords.split("-") #1178914-1248724

    if cluster_name != 'CP026304':
        continue
    
    for cluster_key in id_to_features.keys():
        gbidfull, cluster_num_full = cluster_key.split("|")
        if len(gbidfull.split("_")) > 2:
            gbid = gbidfull.rsplit("_", 2)[0]
            coord_start, coord_end = gbidfull.split("_")[1:] #1054705_1355970
        else:
            gbid = gbidfull
            coord_start, coord_end = (None, None)
        
        if cluster_name == gbid:
            cluster_num = cluster_num_full.split('cluster-')[1]
            if int(cluster_num) != int(out_cluster_num):
                continue

            if coord_start != None:
                if int(coord_start) < int(out_coord_start) and int(out_coord_end) < int(coord_end):
                    print "\nNo coord in gbidfull: ", line.split("\t")[:5], id_to_features[cluster_key]
                    print line,
                    outf.write(line)
                    outf.write("\t")

                    for domain in domains:
                        if domain not in id_to_features[cluster_key].keys():
                            print 0,
                            outf.write("0\t")
                        else:
                            print id_to_features[cluster_key][domain],
                            outf.write(str(id_to_features[cluster_key][domain]))
                            outf.write("\t")
                    outf.write("\n")

            else:
                print "\n Gbid from nt: ", line.split("\t")[:5], id_to_features[cluster_key]
                outf.write(line)
                outf.write("\t")
                print line,

                for domain in domains:
                    if domain not in id_to_features[cluster_key].keys():
                        print 0,
                        outf.write("0\t")
                    else:
                        print id_to_features[cluster_key][domain],
                        outf.write(str(id_to_features[cluster_key][domain]))
                        outf.write("\t")
                outf.write("\n")
                
outf.close()
