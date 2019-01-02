#!/usr/bin/env python
import json
import sys
import glob
import os

# Parse json files and print the genbank ids of only polyketidespython
# parse_json.py > mibig.gbids.pks
# cat mibig.gbids.pks |wc
#     493     493    4581

json_files = glob.glob("mibig_json_1.4/*.json")
for json_file in json_files:
    with open(json_file, 'rU') as f:
        b = json.loads(f.read())
        # print b
        p = b['general_params']
        # print p
        l = p['loci']
        if not 'Polyketide' in p['biosyn_class']:
            continue

        for x in l['nucl_acc']:
            gbid = x['Accession']
            print os.path.basename(json_file).split(".")[0], "\t",  gbid.strip().split('.')[0]
            #print gbid.strip().split('.')[0]

