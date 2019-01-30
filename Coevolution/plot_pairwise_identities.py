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
#
# xyz=np.array(np.random.random((100,3)))
# plt.scatter(xyz[:,0], xyz[:,1])
# plt.savefig('foo.png')
#
# sys.exit(0)


# This code will make a 6deb production plot. The input is 6deb data: PCC<tab>sterror<tab>48htime point)
ks_pair_identities = defaultdict(list)
target_pair_identities = defaultdict(list)

# pairwis_filename = "pairwise_identities.mibig.out"
pairwise_filename = "pairwise_identities.mibig.out.long"
with open(pairwise_filename, 'r') as f:
    data = f.readlines()
    data = map(lambda x: x.strip(), data)
    for line in data:
        print line
        pair, ks_pair_identity, target_pair_identity = line.split("\t")
        # pairwise_identities.12.5kb.out.long
        gene1, gene2 = pair.split("||")
        target = gene1.split("|")[1]
        print gene1, gene2, target
        ks_pair_identities[target].append(float(ks_pair_identity))
        target_pair_identities[target].append(float(target_pair_identity))

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

for target in ks_pair_identities:
    # plt.scatter(ks_pair_identities, target_pair_identities, color='r', s=1)
    if target not in target_to_color:
        continue
    plt.scatter(ks_pair_identities[target], target_pair_identities[target], color=target_to_color[target], s=100)
    plt.xlabel('KS1-KS2 identity', size=10)
    plt.ylabel('Trget1-Target2 identity', size=10)

plt.xlim((0,1))
plt.ylim([0,1])
# plt.savefig('mibig.png')
plt.savefig('12.mibig.long.png', dpi=400)

# fig = plt.figure(5, figsize=(10, 6))
# plt.savefig("image.png",bbox_inches='tight', dpi=100)
#plt.show()
