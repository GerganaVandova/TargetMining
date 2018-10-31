import json
import sys
import glob

# Parse json files and print the genbank ids of only polyketidespython
# parse_json.py > mibig.gbids.pks
# cat mibig.gbids.pks |wc
#     493     493    4581

ids_to_class = {}
ids = []

json_files = glob.glob("mibig_json_1.4/*.json")
for json_file in json_files:
    with open(json_file, 'rU') as f:
        b = json.loads(f.read())
        # print b
        p = b['general_params']
        # print p
        l = p['loci']
        bsc = p['biosyn_class'][0]
        # print bsc
        gbid = l['nucl_acc'][0]['Accession']
        # print gbid
        if bsc == 'Polyketide':
            ids_to_class[gbid] = bsc
            ids.append(gbid)

# for key in ids_to_class.keys():
#     print key

for id in ids:
    print id
