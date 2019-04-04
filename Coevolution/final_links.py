#!/usr/bin/env python
# mv publish/12targets.20kb/ /var/www/html/gvandova
# chmod +x -R f12targets.20kb/

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

def get_antismash_links(summary_filename):

    f = open("Clusters.609.10kb.txt.publish", "w")
    header = open(summary_filename).readlines()[0]
    f.write("%s\n" % header)

    lines = open(summary_filename).readlines()[1:]
    for line in tqdm.tqdm(lines):
        line = line.strip()
        # print line
        params1 = "\t".join(line.split("\t")[:4])
        antismash_link = line.split("\t")[4]
        params2 = "\t".join(line.split("\t")[5:])
        # /Users/gvandova/TargetMining/Antismash_gbids/antismash_output_assemblies_all/CP003233_1726930_2028193/index.html
        if "assemblies" in antismash_link:
            new_link = antismash_link.replace("/Users/gvandova/TargetMining/Antismash_gbids/antismash_output_assemblies_all/", "http://maguro.stanford.edu/gvandova/609targets.10kb/")
        else:
            new_link = antismash_link.replace("/Users/gvandova/TargetMining/Antismash_gbids/antismash_output/", "http://maguro.stanford.edu/gvandova/609targets.10kb/")

        print line, "\n"
        print "%s\t%s\t%s\n" % (params1, new_link, params2)
        f.write("%s\t%s\t%s\n" % (params1, new_link, params2))

    f.close()

def main():

    summary_filename = "Clusters.609.10kb.txt"
    get_antismash_links(summary_filename)



if __name__ == "__main__":
    main()
