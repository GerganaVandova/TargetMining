#!/usr/bin/env python
import subprocess
import os
import os.path
import tqdm
import signal
import sys
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
from multiprocessing import Pool

TEMP_ANTISMASH_OUTPUT = "/home/gvandova/antismash_output_assemblies_2018/"
FINAL_ANTISMASH_OUTPUT = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output_assemblies_part"
GB_FOLDER = "/home/gvandova/assemblies_part"



def f(gbid):
    print "Start %s" % gbid
    temp_dir = os.path.join(TEMP_ANTISMASH_OUTPUT)
    final_dir = os.path.join(FINAL_ANTISMASH_OUTPUT, gbid)

    local_filename = os.path.join(GB_FOLDER, gbid + ".gbff")
    
    geneclusters_filename = os.path.join(final_dir, "geneclusters.js")
    print "checking if %s exists" % geneclusters_filename
    if os.path.exists(geneclusters_filename) == True:
        print geneclusters_filename, " exist"
        return
    
    antismash = "/home/gvandova/bin/run_antismash" + " " + local_filename + " " + temp_dir
    print antismash
    os.system(antismash)
    print "done %s" % gbid
    subprocess.call(["rm", "-rf", final_dir])
    subprocess.call(["mv", "-f", os.path.join(temp_dir, gbid), FINAL_ANTISMASH_OUTPUT])
    print " ".join(["mv", "-f", os.path.join(temp_dir, gbid), FINAL_ANTISMASH_OUTPUT])

def load_gbids():
    gbids = set()
    for f in os.listdir(GB_FOLDER):
        gbid = f.split(".gbff")[0]
        gbids.add(gbid)
    return gbids


if __name__ == "__main__":
    # Get coordinates of KS from Blast results 
    # Run Antismash 
    gbids = load_gbids()
    total = len(gbids)
    print "Total %s" % total
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    p = Pool(processes=45)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        for _ in tqdm.tqdm(p.imap_unordered(f, gbids), total=len(gbids)):
            pass

    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        p.terminate()
    else:
        print("Normal termination")
        p.close()
    p.join()

    print "Done, exiting"

