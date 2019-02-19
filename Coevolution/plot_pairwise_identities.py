#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pyplot
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work
from collections import defaultdict

ks_pair_identities = defaultdict(list)
target_pair_identities = defaultdict(list)
pairs = []

# pairwis_filename = "pairwise_identities.mibig.out"
pairwise_filename = "pairwise_identities.12.5kb.out"
with open(pairwise_filename, 'r') as f:
    data = f.readlines()
    data = map(lambda x: x.strip(), data)
    for line in data:
        print line
        gene1, gene2, len1, len2, ks_pair_identity, target_pair_identity = line.split("\t")
        target = gene1.split("|")[1]
        print gene1, gene2, target
        ks_pair_identities[target].append(float(ks_pair_identity))
        target_pair_identities[target].append(float(target_pair_identity))
        pairs.append((gene1,gene2))

# ks_pair_identities = map(float, ks_pair_identities)
# target_pair_identities = map(float, target_pair_identities)

print len(ks_pair_identities), len(target_pair_identities)

target_to_color = {
         "AdmT_ACC": 'blue',
         "SalI_beta_proteasome": 'lightblue',
         "GriR_DnaN": "cyan",
         "EF-Tu": "midnightblue",
         "PtmP3_FabB-F": "r",
         "BatG_FabI": "darkmagenta",
        "GyrB-R": "magenta",
        "mupM_Ile-tRNA-syn": "brown",
        "borI_Thr-tRNA-syn": "lightpink",
        "agnB2_Leu-tRNA-syn": "orange",
        "rubR1_TIF": "green",
        "Ind0_Trp-tRNA-syn": "lightgreen"}

target_to_name = {
         "AdmT_ACC": 'Acetyl CoA carboxylase',
         "SalI_beta_proteasome": 'Beta proteasome subunit',
         "GriR_DnaN": "DNA polymerase sliding clamp",
         "EF-Tu": "Elongation factor Tu",
         "PtmP3_FabB-F": "3-oxoacyl-[acyl-carrier-protein] synthase 1 (FabB/F)",
         "BatG_FabI": "Enoyl-[acyl-carrier-protein] reductase [NADH] FabI",
        "GyrB-R": "Gyrase B",
        "mupM_Ile-tRNA-syn": "Isoleucyl tRNA synthetase",
        "borI_Thr-tRNA-syn": "Threonyl-tRNA synthetase",
        "agnB2_Leu-tRNA-syn": "Leucil-tRNA synthase",
        "rubR1_TIF": "Translation initiation factor",
        "Ind0_Trp-tRNA-syn": "Tryptophanyl-tRNA synthase"}

for target in ks_pair_identities:
    # plt.scatter(ks_pair_identities, target_pair_identities, color='r', s=1)
    if target not in target_to_color:
        continue
    color = target_to_color[target]
    label = target_to_name[target]
    plt.scatter(ks_pair_identities[target], target_pair_identities[target],
                color=color, s=10, label=label)
    plt.xlabel('KS1-KS2 identity', size=10)
    plt.ylabel('Target1-Target2 identity', size=10)
    plt.legend(loc=4, fontsize="xx-small")


    # y = [2.56422, 3.77284, 3.52623, 3.51468, 3.02199]
    # z = [0.15, 0.3, 0.45, 0.6, 0.75]
    # n = [58, 651, 393, 203, 123]
    #
    # fig, ax = plt.subplots()

# #Doesn't work
# for i, txt in enumerate(pairs):
#     plt.annotate(txt, (ks_pair_identities[i], target_pair_identities[i]))


plt.xlim([0, 100])
plt.ylim([0, 100])
# plt.savefig('mibig.png')
plt.savefig('12.5kb.annotated.png', dpi=400)
