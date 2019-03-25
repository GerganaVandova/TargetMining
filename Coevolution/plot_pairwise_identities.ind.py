#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib import pyplot
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work
from collections import defaultdict
from Bio import SeqIO


def get_targets(antismash_outfilename):
    # Get list of target gene names
    # Head of file:
    # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
    targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target = line.split("|")[1]
        targets.add(target)
    return(targets)


def target_to_name(targets_fasta):
    # Get target name from targets fasta file
    # >DEG10180001_Molybdopterin_biosynthesis_mog_protein
    # >tclQ_L11 tr|A0A097PTA1|A0A097PTA1_9STAP 50S ribosomal protein L11
    target_to_names = defaultdict(list)
    for record in SeqIO.parse(open(targets_fasta, "rU"), "fasta"):
        gbidfull = record.id
        gbid, name = gbidfull.split("_", 1)  # for E. coli 609 targets
        # print gbidfull
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
            print line
            gene1, gene2, len1, len2, ks_pair_identity, target_pair_identity,\
                d = line.split("\t")
            qtarget = gene1.split("|")[1]
            if qtarget != target:
                continue
            # print gene1, gene target

            # Remove identical clusters: ran script for 12 targets with this, when there is only one pair, it wont be displayed
            # TODO Rerun with this comented
            # if float(ks_pair_identity) == 100 and float(target_pair_identity) == 100:
            #     continue
            ks_pair_identities[target].append(float(ks_pair_identity))
            target_pair_identities[target].append(float(target_pair_identity))
    return (ks_pair_identities, target_pair_identities)


def plot_identities(targets, colors, target_to_names, ks_pair_identities, target_pair_identities, qtarget, mode, subplot=None):

    # For 92 nd 609 targets:
    # target_to_color = dict([(k, v) for k, v in zip(targets, colors)])
    # target_to_name = dict([(k, k) for k in targets]) # comment when long target names needed

    # For 12 targets:
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
        print target
        if target == qtarget:
            color = 'black'
            if target not in target_to_color:
                color = "black"
                print "target not in target to color", target
            else:
                color = target_to_color[target]
                print target, target_to_color[target]
            if target not in target_to_names:
                print "target %s not in target to names, label = None" % target
                label = None
            else:
                label = target_to_names[target]
            # print "ks_pair identities", ks_pair_identities, color

            if mode == "subplot":
                color = 'black'
                print subplot
                # (6,6 for 92.5kb, 6,8 for 92.10kb, 8,9 609.5kb, 9,9 609.10kb)
                ax = plt.subplot(4, 3, subplot)
                ax.set_title(label, fontweight="bold", size=5)  # Title size=5 for 36 plots

                # For subplots single target per plot:
                plt.scatter(ks_pair_identities[target],
                            target_pair_identities[target],
                            color=color, s=0.5, label=label)
                plt.xlim([0, 100])
                plt.ylim([0, 100])
                plt.xlabel('KS1-KS2', size=5)
                plt.ylabel('Target1-Target2', size=5)

                # If you don't want scale displayed:
                # ax.set_xticklabels([])
                # ax.set_yticklabels([])

                # If you want scale displayed
                plt.tick_params(axis='x', labelsize=4)  # size=5 for 36 plots
                plt.tick_params(axis='y', labelsize=4)

                # To keep plots square
                ax.set_aspect('equal')

                # plt.subplots_adjust(wspace=0, hspace=0)  # 12.5kb 9 subplots
                plt.subplots_adjust(wspace=0.4, hspace=0.5)  # 12.10kb, 92.5kb 36 subplots
                # plt.subplots_adjust(wspace=0.4, hspace=0.7) # 92.10kb 68 subplots

                # removing the tick marks
                ax.tick_params(top='off', bottom="off", left="off", right='off')

                # Edge of the graph thickness
                ax.spines['left'].set_linewidth(0.5)
                ax.spines['right'].set_linewidth(0.5)
                ax.spines['top'].set_linewidth(0.5)
                ax.spines['bottom'].set_linewidth(0.5)

            if mode == "individual":
                plt.xlabel('KS1-KS2 identity', size=10)
                plt.ylabel('Target1-Target2 identity', size=10)
                plt.xlim([0, 100])
                plt.ylim([0, 100])
                plt.axes().set_aspect('equal')

                # One target per plot:
                # plt.scatter(ks_pair_identities[target],
                #             target_pair_identities[target],
                #             color=color, s=10, label=label)
                # plt.legend(loc=4, fontsize="xx-small")
                # figname = "plots_12/mibig." + qtarget + ".png"
                # plt.savefig(figname, dpi=400)
                # plt.clf()

                # Single plot, all targets:
                plt.scatter(ks_pair_identities[target],
                            target_pair_identities[target],
                            color=color, s=20, label=label) # s = 5 for 12, 92, 609 targets
                # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize="xx-small")
                plt.legend(loc=4, fontsize="xx-small") # for mibig
                figname = "plots_12/mibig." + qtarget + ".all.png"
                plt.savefig(figname, dpi=400, bbox_inches='tight') # so that legend is not cut off when plot square


def main():

    colors = [
        "#000000", "#012C58", "#1CE6FF", "#FF34FF",  # "#FFFF00" yellow, replaced with "#012C58"
        "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC",
        "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#7ED379", "#1B4400",  # "#FEFFE6" yellow, replaced with , "#7ED379"
        "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0",
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

    pairwise_filename = "pairwise_identities.mibig.out"
    antismash_outfilename = "../Antismash_gbids/out.12.filtered.10kb"
    targets_fasta = "../Antismash_gbids/targets.12.fa.cleannames"

    targets = get_targets(antismash_outfilename)
    print "len of targets:", len(targets), "len of colors: ", len(colors)

    target_to_names = target_to_name92(targets_fasta) # for 12 and 92 targets, split in different way
    # target_to_names = target_to_name(targets_fasta) # for 609 targets, split in different way

    # mode = "subplot"
    mode = "individual"
    i = 1
    for target in targets:
        # if target == "DEG10180179": # this is Malonyl_CoA-acyl_carrier_protein_transacylase
        #     continue
        # if target != "DEG10180155":
        #     continue
        ks_pair_identities, target_pair_identities = \
            get_pairs(pairwise_filename, target)

        print len(ks_pair_identities), len(target_pair_identities)
        if len(ks_pair_identities[target]) <= 1:
            continue

        plot_identities(targets,
                        colors,
                        target_to_names,
                        ks_pair_identities,
                        target_pair_identities,
                        target,
                        mode,
                        subplot=i)
        i += 1
        if i > 12:
            break

    if mode != "individual":
        plt.savefig("plots_12/mibig.subplot.png", dpi=400)
        plt.clf()



if __name__ == "__main__":
    main()
