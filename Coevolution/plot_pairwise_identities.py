#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import pyplot
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work

# This code will make a 6deb production plot. The input is 6deb data: PCC<tab>sterror<tab>48htime point)

ks_pair_identites = []
target_pair_identities = []

with open(sys.argv[1], 'r') as f:
    data = f.readlines()
    data = map(lambda x: x.strip(), data)
    for line in data:
        print line
        ks_pair, ks_pair_identity, target_pair, target_pair_identity = line.split("\t")
        ks_pair_identites.append(ks_pair_identity)
        target_pair_identities.append(target_pair_identity)

ks_pair_identites = map(float, ks_pair_identites)
target_pair_identities = map(float, target_pair_identities)

plt.scatter(ks_pair_identites, target_pair_identities)
fig = plt.figure(5, figsize=(10, 6))
plt.savefig("image.png",bbox_inches='tight', dpi=100)
# plt.show()

sys.exit(0)

print pccs, "\n", sterrs, "\n", debs
pos = np.arange(0, 16, 2)    # the bar centers on the y axis for 22 PCCs
print pos
b = plt.barh(pos, debs, xerr=sterrs, align='center', height=1, color = '0.65', ecolor='k', label='48h')
#plt.grid(True)
plt.yticks(pos, pccs)
plt.xlabel('6deb to standard ion counts ratio at 48h', size=15)
#plt.title('6deB production of 22 PCCs')
plt.tick_params(labelsize=14, pad=10)
plt.xlim((0,6)) #for 8 PCCs
#remove pads from both sides of barplot
plt.ylim([-2,16])
plt.tight_layout() #Show long labels
plt.show()
