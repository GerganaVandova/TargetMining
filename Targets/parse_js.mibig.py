#!/usr/bin/env python
from Bio import SeqIO
import json
import os
from pprint import pprint
import glob
import codecs
import sys
import requests
import urllib
import subprocess
UTF8Writer = codecs.getwriter('utf8')
sys.stdout = UTF8Writer(sys.stdout)


def get_targets():
    path = "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Targets/mibig_json.v1.4/*.json"
    # outtargetfile = "mibig_targets.txt"
    fclusters = open("mibig_clusters.txt", "w")
    ftargets = open("mibig_targets.txt", "w")

    targets = []

    files = glob.glob(path)
    for file in files:
        # print file
        f = open(file, 'r')
        data_json = f.read()
        data_js = json.loads(data_json)
        mibigid = data_js["general_params"]["mibig_accession"]
        biosyn_class = data_js["general_params"]["biosyn_class"]
        compounds = data_js["general_params"]["compounds"]
        biosyn_class_clean = (", ".join(biosyn_class))
        for c in compounds:
            compound = c.get("compound")
            target = c.get("chem_target")
            targets.append(target)
            activity = c.get("chem_act", [])
            activities = (", ".join(activity))

            # If you need cluster coordinates:
            # start = data_js["general_params"]["loci"]["nucl_acc"][0]["start_coord"]
            # end = data_js["general_params"]["loci"]["nucl_acc"][0]["end_coord"]

            # print("%s\t%s\t%s\t%s\t%s" %
            # (mibigid, compound, biosyn_class_clean, target, activities))
            fclusters.write("%s\t%s\t%s\t%s\t%s" %
                            (mibigid, compound, biosyn_class_clean,
                             target, activities))
            fclusters.write("\n")
        f.close()

    fclusters.close()

    new_targets = []

    for target in targets:
        try:
            target = target.lower()
            new_targets.append(target)
        except:
            continue

    set_targets = sorted(set(new_targets))
    print(set_targets)
    print(len(set_targets))

    for target in set_targets:
        ftargets.write(str(target))
        ftargets.write("\n")

        # ftargets.write("\n")
    ftargets.close()
    return set_targets


def get_uniprot(target):

    query = target.strip()
    f = {
        'fil': 'reviewed:yes',
        'sort': 'score',
        'columns': 'id,entry name,reviewed,protein names,genes,organism,length',
        'format': 'fasta',
        'query': query
    }
    url = 'https://www.uniprot.org/uniprot/?' + urllib.urlencode(f)
    #url = 'https://www.uniprot.org/uniprot/?fil=reviewed%3Ayes&sort=score&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&format=fasta&query=gyrb%2C+subunit+b+protein+of+dna+gyrase'
    print "target: ", query, "\n", url
    r = requests.get(url, allow_redirects=True)
    query1 = query.replace(",", "")
    query2 = query1.replace("/", "_")
    query_new = query2.replace(" ", "_")
    outfile = "Uniprot/" + query_new + '.uniprot.fasta'
    open(outfile, 'wb').write(r.content)
    records = list(SeqIO.parse(outfile, "fasta"))
    print "number of sequences: ", len(records), "\n\n"

    if len(records) > 500:
        f = {
            'fil': 'identity:0.5',
            'sort': 'score',
            'columns': 'id,entry name,reviewed,protein names,genes,organism,length',
            'format': 'fasta',
            'query': query
        }
        url = 'https://www.uniprot.org/uniref/?' + urllib.urlencode(f)
        print "target: ", query, "\n", url
        r = requests.get(url, allow_redirects=True)
        outfile = "Uniprot/" + query_new + '.uniref.fasta'
        open(outfile, 'wb').write(r.content)
        records = list(SeqIO.parse(outfile, "fasta"))
        print "number of sequences: ", len(records), "\n\n"
        os.remove(os.path.join("Uniprot/", query_new + '.uniprot.fasta'))


targets = get_targets()
for target in targets:
    if target == "":
        continue
    # target = 'gyrb, subunit b protein of dna gyrase'
    get_uniprot(target)
