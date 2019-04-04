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
import tqdm
import json
import subprocess

def copy_antismash(summary_filename):

    lines = open(summary_filename).readlines()[1:]
    for line in tqdm.tqdm(lines):
        line = line.strip()
        # print line
        antismash_link = line.split("\t")[4]
        # /Users/gvandova/TargetMining/Antismash_gbids/antismash_output_assemblies_all/CP003233_1726930_2028193/index.html
        new_link = antismash_link.replace("/Users/gvandova/", "/mnt/gnpn/gnpn/projects/orphanpks/")
        # print new_link
        path = new_link.split("/index.html")[0]
        print path
        dest_folder = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/publish/609targets.10kb"
        is_assemblydir = path.split("/")[8]
        # if "assemblies_all" in is_assemblydir:
        #     dest_folder = "/home/gvandova/publish/antismash_output_assemblies"
        # else:
        #     dest_folder = "/home/gvandova/publish/antismash_output_nt"
        # print is_assemblydir
        # print path
        subprocess.call(['cp', '-a', path, dest_folder])


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

    # summary_filename = "Clusters.609.10kb.txt"
    # copy_antismash(summary_filename)

    summary_filename = "Clusters.12.20kb.txt"
    copy_genbank_files(summary_filename)



if __name__ == "__main__":
    main()
