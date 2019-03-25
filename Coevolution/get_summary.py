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


def get_targets(antismash_outfilename):
    # Get list of target gene names
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
    target_to_names = defaultdict(list)
    for record in SeqIO.parse(open(targets_fasta, "rU"), "fasta"):
        gbidfull = record.id
        gbid, name = gbidfull.split("_", 1)  # for E. coli 609 targets
        name = name.replace("/", " ")
        name = name.replace("_", " ")
        target_to_names[gbid] = name
    return target_to_names


def target_to_name92(targets_fasta):
    # Get target name from targets fasta file of 12 targets and 92 targets lists
    # >tclQ_L11 tr|A0A097PTA1|A0A097PTA1_9STAP 50S ribosomal protein L11
    target_to_names92 = defaultdict(list)
    targets_file = open(targets_fasta).readlines()
    for line in targets_file:
        if not line.startswith(">"):
            continue
        line = line.strip()
        target_id, name = line.split(" ", 1)  # for 92 targets
        target_id = target_id.split(">")[1]
        # name = name.replace("/", " ")
        name = name.replace("_", " ")
        target_to_names92[target_id] = name
    return target_to_names92


def get_species(speciesfilename):
    # Get species name from species file
    # BDBI01000023   Nocardia sp.
    gbid_to_species = defaultdict()

    speciesfile = open(speciesfilename).readlines()
    for line in speciesfile:
        gbid, name = line.split("\t")
        gbid_to_species[gbid] = name
    return gbid_to_species


def get_phyla(taxafilename):
    # Get phyla from taxa file
    # BDBI01000023   Bacteria    Actinobacteria  Corynebacteriales   Nocardiaceae    Nocardia
    gbid_to_phyla = defaultdict()
    taxafile = open(taxafilename).readlines()
    for line in taxafile:
        feats = line.split("\t", 3)
        if len(feats) < 3:
            continue
        gbid, taxa, phyla, rest = feats
        gbid_to_phyla[gbid] = phyla
    return gbid_to_phyla


def get_coevolution_score(coevolutionfilename):
    gbid_to_score = defaultdict()
    coevolutionfile = open(coevolutionfilename).readlines()
    for line in coevolutionfile:
        gbid, target_id, count, score = line.split("\t")
        mytuple = (target_id, count, score)
        value = "|".join(mytuple)
        gbid_to_score[gbid] = value
    return gbid_to_score


def get_antismash_link(gbid, ks_start, ks_end):
    gbid_to_antismashlink = defaultdict()
    # Local path
    gbdir_nt = "/Users/gvandova/TargetMining/Antismash_gbids/antismash_output"
    gbdir_as = "/Users/gvandova/TargetMining/Antismash_gbids/antismash_output_assemblies_all"
    path_nt = os.path.join(gbdir_nt, gbid, "index.html")
    path_as = os.path.join(gbdir_as, gbid, "index.html")

    # Maguro path
    gbdir1 = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output"
    gbdir2 = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output_assemblies_all"
    path1 = os.path.join(gbdir1, gbid, "index.html")
    path2 = os.path.join(gbdir2, gbid, "index.html")

    if os.path.exists(path1) is True:
        antismashfile = path_nt

    elif os.path.exists(path2) is True:
        antismashfile = path_as

    else:
        folders = []
        for i in os.listdir(gbdir2):
            if os.path.exists(os.path.join(gbdir2, i)) and gbid in i:
                folders.append(i)
        for folder in folders:
            gbid_folder, gbid_start, gbid_end = folder.split("_")
            gbid_start = int(gbid_start)
            gbid_end = int(gbid_end)
            if gbid_start < ks_start and gbid_end > ks_end:
                antismashfile = os.path.join(gbdir_as, folder, "index.html")
    gbid_to_antismashlink[gbid] = antismashfile
    return gbid_to_antismashlink


def get_second_copy(second_copy_filename):
    # agnB2_Leu-tRNA-syn	JOID01000025	2	136910	False	67565	70387	307	813	137570	1e-126
    gbid_to_copy = defaultdict()
    copyfile = open(second_copy_filename).readlines()
    for line in copyfile:
        target_id, gbid, copynum, genome_size, complete, rest = line.split("\t", 5)
        mytuple = (copynum, genome_size, complete)
        value = "|".join(mytuple)
        gbid_to_copy[gbid] = value
    return gbid_to_copy


def main():
    # targets_fasta = "../Antismash_gbids/targets.92.fa.cleannames"
    # target_to_names = target_to_name92(targets_fasta)

    targets_fasta = "../Antismash_gbids/targets.609.fa.longnames"
    target_to_names = target_to_name(targets_fasta)

    speciesfilename = "../Genbank/species.txt"
    gbid_to_species = get_species(speciesfilename)

    taxafilename = "../Genbank/taxa.txt"
    gbid_to_phyla = get_phyla(taxafilename)

    coevolutionfilename = "coevolution_scores.609.10kb"
    gbid_to_score = get_coevolution_score(coevolutionfilename)

    second_copy_filename = "../Second_copy/out.second_copy.609.10kb.filtered"
    gbid_to_copy = get_second_copy(second_copy_filename)

    f = open("Clusters.609.10kb.txt", "w")
    f.write("target_name\tgenbankid\ttargetid\tks_start\tks_end\target_start\ttarget_end\tcluster_num\tcluster_type\tcluster_coord\tcluster_len\tspecies\tphyla\tcount\tscore\tcopynum\tgenome_size\tcomplete\tantismashfile\n")
    # Read KS fasta file and write phyla for each gbid

    for record in tqdm.tqdm(SeqIO.parse(open("KS.609.10kb.fasta", "rU"), "fasta")):
        gbidfull = record.id
        seq = record.seq
        gbid, target_id, \
            ks_start, ks_end, \
            target_start, target_end, \
            cluster_num, cluster_type, \
            cluster_coord, tandem_dist = gbidfull.split("|")
        ks_start = int(ks_start)
        ks_end = int(ks_end)

        target_name = target_to_names[target_id]
        descr = "\t".join(gbidfull.split("|"))

        gbid_to_antismashlink = get_antismash_link(gbid, ks_start, ks_end)
        antismashfile = gbid_to_antismashlink[gbid]

        if gbid not in gbid_to_phyla.keys():
            phyla = "None"
        else:
            phyla = gbid_to_phyla[gbid]

        if gbid not in gbid_to_species.keys():
            species = "None"
        else:
            species = gbid_to_species[gbid]
            species = species.strip()

        if gbid not in gbid_to_score.keys():
            coevolution_score = target_id + '|0|' + '-1'
        else:
            coevolution_score = gbid_to_score[gbid]
        target_ids, count, score = coevolution_score.split("|")
        score = float(score)

        second_copy = gbid_to_copy[gbid]
        copynum, genome_size, complete = second_copy.split("|")

        # print "%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s" % \
        #     (target_name, descr, species, phyla, count, score,
        #      copynum, genome_size, complete, antismashfile)
        f.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\n" %
                (target_name, descr, species, phyla, count, score,
                 copynum, genome_size, complete, antismashfile))
    f.close()

if __name__ == "__main__":
    main()
