#!/usr/bin/env python
import matplotlib
import numpy as np
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
# import pandas as pd
from matplotlib import pyplot
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work
from collections import defaultdict
from Bio import SeqIO
import json
import tqdm
import subprocess

def copy_antismash(summary_filename):
    lines = open(summary_filename).readlines()[1:]
    for line in tqdm.tqdm(lines):
        line = line.strip()
        antismash_link = line.split("\t")[12]
        print antismash_link
        # /Users/gvandova/TargetMining/Antismash_gbids/antismash_output_assemblies_all/CP003233_1726930_2028193/index.html
        new_link = antismash_link.replace("/Users/gvandova/",
                                          "/mnt/gnpn/gnpn/projects/orphanpks/")


        print new_link
        path = new_link.split("/index.html")[0]
        # full_path = 'maguro:/' + path
        print path
        dest_folder = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/publish/616targets.10kb"
        # dest_folder = "/Users/gvandova/Dropbox/publish/12targets/"
        subprocess.call(['cp', '-a', path, dest_folder])
        # subprocess.call(['scp', '-r', full_path, dest_folder])


def copy_genbank_files(summary_filename):
    f = open(summary_filename).readlines()[1:]
    for line in f:

        line = line.strip()
        gbid = line.split("\t")[2]
        # if "CP011497" not in gbid:
        #     continue
        clusternum = line.split("\t")[5]  #cluster-1
        n = clusternum.split("cluster-")[1]
        antismash_link = line.split("\t")[4]
        new_link = antismash_link.replace("/Users/gvandova/", "/mnt/gnpn/gnpn/projects/orphanpks/")
        # print new_link
        path = new_link.split("/index.html")[0]
        print path
        for file in os.listdir(path):
            if str(n)+".gbk" in file:
                print file
                gb_filename = file
        full_gb_filename = os.path.join(path, gb_filename)
        print full_gb_filename
        dest_folder = "/home/gvandova/bin/example/antismash_gbks"
        subprocess.call(['cp', '-a', full_gb_filename, dest_folder])

def main():

    summary_filename = "Clusters.616.10kb.txt"
    print summary_filename
    copy_antismash(summary_filename)

    # summary_filename = "Clusters.12.10kb.txt"
    # copy_genbank_files(summary_filename)



if __name__ == "__main__":
    main()
