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
import scipy.stats
import seaborn as sns
from textwrap import wrap


# Calculating p-value to display on graphs
def p_value_def(p_value):

    if p_value < 0.001:
        return_p_val = r"$p<0.001$"
    elif p_value < 0.01:
        return_p_val = r"$p<0.01$"
    elif p_value < 0.05:
        return_p_val = r"$p<0.05$"
    else:
        return_p_val = r"ns"
    return return_p_val

pairwise_filenames = ["pairwise_identities.mibig.out",
                      "pairwise_identities.14.5kb.out",
                      "pairwise_identities.14.10kb.out",
                      "pairwise_identities.14.m50kb.out"]
i = 1
s = ['a) ', 'b) ', 'c) ', 'd) ', 'e) ', 'f) ']
for pairwise_filename in pairwise_filenames:
    try:
        name = pairwise_filename.split('pairwise_identities.14.')[1]
        distance = name.split('.out')[0]
        title = s[i-1] + "14 targets " + distance + " cutoff"
    except:
        # name = pairwise_filename.split('mibig.')[1]
        title = s[i-1] + "MIBiG positive set"
    ks_pair_identities = defaultdict(list)
    target_pair_identities = defaultdict(list)
    with open(pairwise_filename, 'r') as f:
        data = f.readlines()
        data = map(lambda x: x.strip(), data)
        for line in data:
            gene1, gene2, len1, len2, ks_pair_identity, target_pair_identity, d = line.split("\t")
            d = float(d)
            target = gene1.split("|")[1]
            ks_pair_identities[target].append(float(ks_pair_identity))
            target_pair_identities[target].append(float(target_pair_identity))
    x = []
    y = []

    color = "black"
    for target in ks_pair_identities:
        x.extend(ks_pair_identities[target])
        y.extend(target_pair_identities[target])

    # Setting a style of figure and font
    sns.set(style="whitegrid", font_scale=1)
    font = {'family': 'sans-serif', 'color': 'black',
            'weight': 'normal', 'size': 8, 'rotation': 0,
            'verticalalignment': 'bottom', 'horizontalalignment': 'left'}

    # calculating linear regression
    # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    # R_square = np.round(r_value * r_value, 2)

    # calculating linear regression with Spearman correlation coefficient
    sp_corr, p_value = scipy.stats.spearmanr(x, y)
    Sp_corr = np.round(sp_corr, 2)
    p_val = np.round(p_value, 2)
    print Sp_corr, p_val

    ax = plt.subplot(2, 2, i)
    plt.subplots_adjust(wspace=0.4, hspace=0.5)  # 12.10kb, 92.5kb 36 subplots
    ax.set_aspect('equal')  # To keep plots square

    # Plot params
    ax.set_title(title, fontweight="bold", size=10)
    plt.grid(False)
    plt.xlabel('KS1-KS2', size=10)  # X Labels
    plt.ylabel('Protein1-Protein2', size=10)  # Y labels
    plt.xlim([0, 100])  # X range
    plt.ylim([0, 100])  # Y range
    plt.tick_params(axis='x', labelsize=10)  # scale size
    plt.tick_params(axis='y', labelsize=10)
    plt.tick_params(top='off', bottom="off",
                    left="off", right='off')  # removing tick marks

    # Plot linear regression
    sns.regplot(x, y,
                fit_reg=False,
                color='grey',
                ci=None,
                line_kws={'linewidth': 1},
                scatter_kws={'s': 3},)

    # Show Rsquare and p-value on plot
    # ax.text(68, 3, r'$R^{2}=$'+'{}'.format(R_square)+"\n" +
    #         p_value_def(p_value), fontdict=font)

    ax.text(68, 3, 'r ='+ '{}'.format(Sp_corr) + "\n" +
            p_value_def(p_value), fontdict=font)

    i += 1
figname = "plots/mibig.correlation.png"
plt.savefig(figname, dpi=400)
