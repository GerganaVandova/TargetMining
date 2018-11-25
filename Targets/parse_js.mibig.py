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
import time
UTF8Writer = codecs.getwriter('utf8')
sys.stdout = UTF8Writer(sys.stdout)


def get_targets():
    path = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Targets/mibig_json.v1.4/*.json"
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
    #print(set_targets)
    #print(len(set_targets))

    for target in set_targets:
        ftargets.write(str(target))
        ftargets.write("\n")

        # ftargets.write("\n")
    ftargets.close()
    return set_targets


def get_uniprot(target):
    query = target.strip()
    if query in ["not determined"]:
        return
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
    query1 = query.replace(",", "")
    query2 = query1.replace("/", "_")
    query3 = query2.replace("(", "_")
    query4 = query3.replace(")", "_")
    query_new = query4.replace(" ", "_")
    outfile_uniprot = "Uniprot/" + query_new + '.uniprot.fasta'
    
    if os.path.exists(outfile_uniprot) == True:
        print outfile_uniprot, " exist\n"
        return
    
    outfile_uniref = "Uniprot/" + query_new + '.uniref.fasta'
    if os.path.exists(outfile_uniref) == True:
        print outfile_uniref, " exist\n"
        return
    
    r = requests.get(url, allow_redirects=True)
    open(outfile_uniprot, 'wb').write(r.content)
    records = list(SeqIO.parse(outfile_uniprot, "fasta"))
    print "number of sequences: ", len(records), "\n\n"
    
    #if len(records) == 0:
    #    os.remove(os.path.join("Uniprot/", query_new + '.uniprot.fasta'))

    if len(records) > 500:
        query_reviewed = query + '+AND reviewed:yes'
        f = {
            'fil': 'identity:0.5',
            'sort': 'score',
            'columns': 'id,entry name,reviewed,protein names,genes,organism,length',
            'format': 'fasta',
            'query': query_reviewed
        }
        
        url = 'https://www.uniprot.org/uniref/?' + urllib.urlencode(f)
        print "target: ", query, "\n", url
        r = requests.get(url, allow_redirects=True)
        open(outfile_uniref, 'wb').write(r.content)
        records = list(SeqIO.parse(outfile_uniref, "fasta"))
        print "number of sequences: ", len(records), "\n\n"
        os.remove(os.path.join("Uniprot/", query_new + '.uniprot.fasta'))
    time.sleep(1)     

#target = 'actin'
#get_uniprot(target)

targets = get_targets()
for target in targets:
    if target == "":
        continue
#    target = '50s_ribosomal_subunit'
    get_uniprot(target)
