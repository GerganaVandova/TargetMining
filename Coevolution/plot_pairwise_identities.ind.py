#!/usr/bin/env python
import matplotlib
import numpy as np
import seaborn as sns
from scipy.stats import spearmanr
import scipy.stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import pyplot
import pandas as pd
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work
from collections import defaultdict
from Bio import SeqIO
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


def get_targets(antismash_ksfasta):
    targets = set()
    for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
        gbidfull = record.id
        target = gbidfull.split("|")[1]
        targets.add(target)
    return(targets)


def target_to_name(targets_fasta):
    # Get target name from targets fasta file
    # >DEG10180001_Molybdopterin_biosynthesis_mog_protein
    # >tclQ_L11 tr|A0A097PTA1|A0A097PTA1_9STAP 50S ribosomal protein L11
    target_to_names = defaultdict(list)
    for record in SeqIO.parse(open(targets_fasta, "rU"), "fasta"):
        gbidfull = record.id
        print gbidfull
        gbid, name = gbidfull.split("_", 1)  # for E. coli 609 targets
        name = name.replace("/", " ")
        name = name.replace("_", " ")
        # print gbid, name
        target_to_names[gbid] = name
    return target_to_names


def target_to_name92(targets_fasta):
    # Get target name from targets fasta file
    # >tclQ_L11 tr|A0A097PTA1|A0A097PTA1_9STAP 50S ribosomal protein L11
    # >DEG10180001_Molybdopterin_biosynthesis_mog_protein
    target_to_names92 = defaultdict(list)
    targets_file = open(targets_fasta).readlines()
    for line in targets_file:
        if not line.startswith(">"):
            continue
        # print line
        line = line.strip()
        gbid, name = line.split(" ", 1)  # for 92 targets
        gbid = gbid.split(">")[1]
        # name = name.replace("/", " ")
        name = name.replace("_", " ")
        # print name
        target_to_names92[gbid] = name
    return target_to_names92


def get_pairs(pairwise_filename, target):

    ks_pair_identities = defaultdict(list)
    target_pair_identities = defaultdict(list)
    with open(pairwise_filename, 'r') as f:
        data = f.readlines()
        data = map(lambda x: x.strip(), data)
        for line in data:
            # print line
            gene1, gene2, len1, len2, ks_pair_identity, target_pair_identity,\
                d = line.split("\t")
            qtarget = gene1.split("|")[1]
            if float(ks_pair_identity) == 100 and float(target_pair_identity) == 100:
                continue
            if qtarget != target:
                continue
            # print gene1, gene2, target
            ks_pair_identities[target].append(float(ks_pair_identity))
            target_pair_identities[target].append(float(target_pair_identity))
    return (ks_pair_identities, target_pair_identities)


def plot_identities(targets, colors, target_to_names, ks_pair_identities,
                    target_pair_identities, qtarget, color, mode, corr_file, subplot=None):

    # For 12 targets
    target_to_color = {
        "AdmT_ACC": "blue",
        "SalI_beta_proteasome": 'lightblue',
        "GriR_DnaN": "cyan",
        "EF-Tu": "midnightblue",
        "PtmP3_FabB-F": "r",
        "BatG_FabI": "darkmagenta",
        "GyrB-R": "magenta",
        "mupM_Ile-tRNA-syn": "brown",
        "borI_Thr-tRNA-syn": "lightpink",
        "agnB2_Leu-tRNA-syn": "orange",
        "methionine_aminopeptidase": "green",
        "Ind0_Trp-tRNA-syn": "lightgreen",
        "ArgK_OTCase": "deeppink",
        "seryl-tRNA_synthetase": "slateblue" }

    target_to_name = {
        "AdmT_ACC": "Acetyl-CoA carboxyltransferase beta subunit",
        "SalI_beta_proteasome": 'Bet proteasome subunit',
        "GriR_DnaN": "DNA polymerase sliding clamp",
        "EF-Tu": "Elongation factor Tu",
        "PtmP3_FabB-F": "3-oxoacyl-[acyl-carrier-protein] synthase 1",
        "BatG_FabI": "Enoyl-[acyl-carrier-protein] reductase [NADH]",
        "GyrB-R": "Gyrase B subunit",
        "mupM_Ile-tRNA-syn": "Isoleucyl-tRNA synthetase",
        "borI_Thr-tRNA-syn": "Threonyl-tRNA synthase",
        "agnB2_Leu-tRNA-syn": "Leucyl-tRNA synthase",
        "methionine_aminopeptidase": "Methionine aminopeptidase",
        "ArgK_OTCase": "Ornithine carbamoyl transferase",
        "seryl-tRNA_synthetase": "Seryl-tRNA synthetase",
        "Ind0_Trp-tRNA-syn": "Tryptophanyl-tRNA synthase"}


    # For 92 and 609 targets:
    target_to_color = dict([(k, v) for k, v in zip(targets, colors)])
    #target_to_color = dict([(k, v) for k, v in zip(ks_pair_identities, colors)])
    target_to_name = dict([(k, k) for k in targets])

    # print "ks_pair_identites", len(ks_pair_identities), ks_pair_identities
    # print "targets", len(targets), targets
    # sys.exit(0)
    # Setting a style of figure and font
    # size = 5 for subplots and 10 for indiv plots
    sns.set(style="whitegrid", font_scale=10)
    font = {'family': 'sans-serif', 'color': 'black',
            'weight': 'normal', 'size': 10, 'rotation': 0,
            'verticalalignment': 'bottom', 'horizontalalignment': 'left'}

    #i = 1
    #for target in target_to_color.keys():
    #    print "a%s<-desc.reordered=='%s'" % (i, target)
    #     i += 1
    # i = 1
    # for target in target_to_color.keys():
    #     print "myCols[a%s]='%s'" % (i, target_to_color[target])
    #     i += 1
    # sys.exit(0)

    n = 0
    for target in ks_pair_identities:
        if target == qtarget:
        # """
        #     if target not in target_to_color:
        #         color = "black"
        #         print "target not in target to color", target, n
        #         n += 1
        #     else:
        #         color = target_to_color[target]
        #         """
            if target not in target_to_names:
                print "target %s not in target to names, label = None" % target
                label = None
            else:
                label = target_to_names[target]

            x = ks_pair_identities[target]
            y = target_pair_identities[target]

            # calculating linear regression
            slope, intercept, r_value, p_value1, std_err = scipy.stats.linregress(x, y)
            R_square = np.round(r_value * r_value, 2)
            p_val1 = np.round(p_value1, 2)

            sp_corr, p_value2 = scipy.stats.spearmanr(x, y)
            Sp_corr = np.round(sp_corr, 2)
            p_val2 = np.round(p_value2, 2)

            corr_file.write("%s\t%s\t%s\t%s\t%s\n" % (target,
                                                          R_square, p_val1,
                                                          Sp_corr, p_val2))

            if mode == "subplot":
                print target, R_square, p_val1, Sp_corr, p_val2
                ax = plt.subplot(3, 4, subplot)
                plt.subplots_adjust(wspace=0.4, hspace=0.5)  # 12.10kb, 92.5kb 36 subplots

                # Plot params
                ax.set_title("\n".join(wrap(label, 30)), fontweight="bold", size=5)
                for i in ['left', 'right', 'top', 'bottom']:  # Edge thickness
                    ax.spines[i].set_linewidth(0.5)

                plt.setp(ax.lines, linewidth=0.1)
                plt.grid(False)
                plt.xlabel('KS1-KS2', size=5)  # X Labels
                plt.ylabel('Protein1-Protein2', size=5)  # Y labels
                plt.xlim([0, 100])  # X range
                plt.ylim([0, 100])  # Y range
                plt.tick_params(axis='x', labelsize=5)  # scale size
                plt.tick_params(axis='y', labelsize=5)
                plt.tick_params(top='off', bottom="off",
                               left="off", right='off')  # removing tick marks

                # Plotting linear regression
                if len(x) == 1:
                    sns.regplot(x, y, fit_reg=False,
                                color='grey',
                                line_kws={'linewidth': 1},
                                scatter_kws={'s': 3},  # 3 for subplot, 10 for indiv plot
                                ci=None)
                else:
                    sns.regplot(x, y, fit_reg=False,  # don't plot linear regression if calculating Spearman correlation
                                color='grey',
                                line_kws={'linewidth': 1},
                                scatter_kws={'s': 3},  # 3 for subplot, 10 for indiv plot
                                ci=None)

                    # ax.text(5, 80, r'$R^{2}=$'+'{}'.format(R_square) + "\n" +
                    #         p_value_def(p_value1),
                    #         'r='+ '{}'.format(Sp_corr) + "\n" +
                    #         p_value_def(p_value2), fontdict=font)
                # Print only Spearman
                ax.text(5, 80, '$r=$'+'{}'.format(Sp_corr) + "\n" +
                        p_value_def(p_value2), fontdict=font)


            if mode == "individual":

                plt.figure(1)
                fig, ax = plt.subplots(figsize=(5, 5))

                # Plot params
                ax.set_title("\n".join(wrap(label, 60)), fontweight="bold", size=10)
                plt.grid(False)
                plt.xlabel('KS1-KS2', size=10)  # X Labels
                plt.ylabel('Protein1-Protein2', size=10)  # Y labels
                plt.xlim([0, 100])  # X range
                plt.ylim([0, 100])  # Y range
                plt.tick_params(axis='x', labelsize=10)  # scale size
                plt.tick_params(axis='y', labelsize=10)
                plt.tick_params(top='off', bottom="off",
                                left="off", right='off')  # removing tick marks

                # color = 'grey'  # for 609 targets

                # Plot linear regression
                if len(x) == 1:
                    sns.regplot(x, y, fit_reg=False, color=color)
                else:
                    sns.regplot(x, y, fit_reg=False, color=color, ci=None)

                    # Print both Pearson and Spearman
                    # ax.text(5, 80, '$Pearson=$'+'{}'.format(R_square) + "\n" +
                    #         p_value_def(p_value1) + "\n" +
                    #         '$Spearman=$'+'{}'.format(Sp_corr) + "\n" +
                    #         p_value_def(p_value2), fontdict=font)

                    # Print only Spearman
                    ax.text(5, 80, '$r=$'+'{}'.format(Sp_corr) + "\n" +
                            p_value_def(p_value2), fontdict=font)

                figname = "plots/616." + qtarget + ".png"
                plt.savefig(figname, dpi=400)
                plt.clf()


def main():

    ci = 0
    colors = [
        "#000000", "#012C58", "#1CE6FF", "#FF34FF",  # "#FFFF00" yellow, replaced with "#012C58"
        "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC",
        "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#7ED379", "#1B4400",  # "#FEFFE6" yellow, replaced with , "#7ED379"
        "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#11915A", "#BA0900", "#6B7900", "#00C2A0",
        "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
        "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF",
        "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF",
        "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
        "#BEC459", "#456648", "#0086ED", "#886F4C",

        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C",
        "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF",
        "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837",
        "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F",
        "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
        "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
        "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C",
        "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3",
        "#CFCDAC", "#D0AC94"        ## added above "#7ED379", "#012C58"
    ]

    pairwise_filename = "pairwise_identities.616.10kb.out"
    antismash_ksfasta = "../Antismash_gbids/KS.616.10kb.fasta.cdhit.90"
    targets_fasta = "../Antismash_gbids/targets.616.fa.cleannames"

    targets = get_targets(antismash_ksfasta)
    print "len of targets:", len(targets), "len of colors: ", len(colors)

    # target_to_names = target_to_name92(targets_fasta)  # for 12 and 92 targets, split in different way
    target_to_names = target_to_name(targets_fasta)  # for 609 targets, split in different way

    # mode = "subplot"
    mode = "individual"
    corr_file = open("616.correlations", 'w')
    i = 0
    for target in targets:
        # if target != "AdmT_ACC":
        #     continue
        print "\n\n**********", target
        ks_pair_identities, target_pair_identities = \
            get_pairs(pairwise_filename, target)

        print len(ks_pair_identities), len(target_pair_identities)
        # print ks_pair_identities, target_pair_identities
        if len(ks_pair_identities[target]) < 1:
            continue
        c = colors[ci]
        ci = (ci + 1) % len(colors)
        plot_identities(targets,
                        colors,
                        target_to_names,
                        ks_pair_identities,
                        target_pair_identities,
                        target,
                        c,
                        mode,
                        corr_file,
                        subplot=i % 12 + 1)

        i += 1
        if i % 12 == 0:
            if mode != "individual":
                plt.savefig("plots/616.10kb.subplots.correlation_part%s.png" % i, dpi=400)
                plt.clf()  # comment if only 12 plots

    if mode != "individual":
        plt.savefig("plots/616.10kb.subplots.correlation_part%s.png" % i, dpi=400)
        plt.clf()

if __name__ == "__main__":
    main()
