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

def get_target_counts(antismash_ksfasta):
    target_counts = {}
    targets = set()
    for record in SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta"):
        gbidfull = record.id
        target = gbidfull.split("|")[1]
        if target in targets:
            target_counts[target] += 1
        else:
            targets.add(target)
            target_counts[target] = 1
    return (targets, target_counts)


def get_targets(antismash_outfilename):
    # Get list of target gene names
    # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|616512-116691|5386
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


def get_coevolution_score_and_homologs(coevolutionfilename):
    gbid_to_score = defaultdict()
    gbid_to_homologs = defaultdict()
    coevolutionfile = open(coevolutionfilename).readlines()
    for line in coevolutionfile:
        # UHIQ01000001	EF-Tu	16	24.5248229672	0	[]	[]	0.18	0.0	0.24	0.01
        gbid, target_id, count, score, homologs, homologs_list, pident_list, R_square, p_val1, Sp_corr, p_val2, \
            = line.strip().split("\t")
        # print homologs_list, len(homologs_list)
        mytuple = (count, Sp_corr, p_val2, score)
        pident_list = pident_list.strip()
        value = "|".join(mytuple)
        gbid_to_score[(gbid, target_id)] = value
        mytuple2 = (homologs, homologs_list, pident_list)
        value2 = "|".join(mytuple2)
        gbid_to_homologs[(gbid, target_id)] = value2

    return (gbid_to_score, gbid_to_homologs)


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


def get_gene_name(genes_fasta):
    # Get target name from targets fasta  of 616 targets file
    # >AM889285|cluster-1|2454838-2501896|t1pks|2500895-2502196|ctg1_139
    gbid_to_gene = defaultdict()
    for record in SeqIO.parse(open(genes_fasta, "rU"), "fasta"):
        gbidfull = record.id
        gbid, _, _, _, gene_coord, gene_name = gbidfull.split("|")
        gene_start, gene_end = gene_coord.split("-")
        gbid_to_gene[gbid, gene_start, gene_end] = gene_name
    return gbid_to_gene


def main():
    targets_fasta = "../Antismash_gbids/targets.119.fa.cleannames"
    target_to_names = target_to_name92(targets_fasta)
    # target_to_names = target_to_name(targets_fasta)

    antismash_ksfasta = "../Antismash_gbids/KS.119.10kb.fasta.cdhit.90"
    targets, target_counts = get_target_counts(antismash_ksfasta)

    speciesfilename = "../Genbank/species.txt"
    gbid_to_species = get_species(speciesfilename)

    taxafilename = "../Genbank/taxa.txt"
    gbid_to_phyla = get_phyla(taxafilename)

    coevolutionfilename = "119.10kb.scores"
    gbid_to_score, gbid_to_homologs = get_coevolution_score_and_homologs(coevolutionfilename)

    second_copy_filename = "../Second_copy/out.second_copy.119.10kb.filtered"
    gbid_to_copy = get_second_copy(second_copy_filename)

    genes_fasta = "../Antismash_gbids/cluster_genes.all.fasta"
    gbid_to_gene = get_gene_name(genes_fasta)

    f = open("Clusters.119.10kb.txt", "w")

    f.write("target_name\ttarget id\tgbid\tscore_final\t")
    f.write("distance\thomologs\thomologs list\ttarget ubiquity\tcopy number\tSpearman correlation\tp value\tcoevolution coefficient\t")
    f.write("antismashfile\tcluster_num\tcluster_type\tgene_name\tspecies\tphyla\ttrue positive\ttrue positives homologs\t")
    f.write("score distance\tscore homologs\tscore target ubiquity\tscore copy number\tscore coevolution\n")

    # Read KS fasta file and write phyla for each gbid
    for record in tqdm.tqdm(SeqIO.parse(open(antismash_ksfasta, "rU"), "fasta")):
        gbidfull = record.id
        seq = record.seq
        gbid, target_id, \
            ks_start, ks_end, \
            target_start, target_end, \
            cluster_num, cluster_type, \
            cluster_coord, d = gbidfull.split("|")
        ks_start = int(ks_start)
        ks_end = int(ks_end)
        d = int(d)
        #
        # if target_id != "AdmT_ACC":
        #     continue
        # if gbid != "FQ859185":
        #     continue

        gene_name = gbid_to_gene[gbid, target_start, target_end]
        target_name = target_to_names[target_id]
        # >CSTD01000001|mupM_Ile-tRNA-syn|1364568|1364989|1356957|1360097|cluster-3|t1pks-nrps|1344560-1399261|4471
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

        coevolution_score = gbid_to_score[(gbid, target_id)]
        count, Sp_corr, p_val2, score = coevolution_score.split("|")
        Sp_corr = float(Sp_corr)
        p_val2 = float(p_val2)

        target_count = target_counts[target_id]
        score = float(score)
        count = int(count)

        # if (gbid, target_id) not in gbid_to_homologs.keys():
        #     print "there are no homologs from the scores file"
        #     homologs = None
        #     homologs_list = None
        #     pident_list = None
        # else:
        homologs, homologs_list, pident_list = gbid_to_homologs[(gbid, target_id)].split("|")
        homologs = int(homologs)
        homologs_list = json.loads(homologs_list)

        second_copy = gbid_to_copy[gbid]
        copynum, genome_size, complete = second_copy.split("|")
        copynum = int(copynum)

        # 1. Calclulate KS-target distance score
        score_distance = -1000
        if d <= 5000:
            score_distance = 1
        else:
            score_distance = 0

        # 2. Calculate score based on occurences of targets
        score_target_frequency = -1000
        if count > 1 and count <= 20:
            score_target_frequency = 1
        else:
            score_target_frequency = 0

        # 3. Calclulate score based on KS-target coevolution
        score_coevolution = -1000
        if Sp_corr > 0.5 and p_val2 < 0.05:
            score_coevolution = 1
        else:
            score_coevolution = 0

        # 4. Calclulate score based on second copy of target
        score_copy = -1000
        if copynum > 1:
            score_copy = 1
        else:
            score_copy = 0

        # 5. Calclulate score based on number of homologs >80% but <90% identity
        score_homologs = -1000
        if homologs >= 1:
            score_homologs = 1
        else:
            score_homologs = 0

        score_final = score_distance + \
            score_target_frequency + \
            score_coevolution + \
            score_copy + \
            score_homologs

        pos_set = set(["KT362046", "LN879418", "JPRX01000001", "KE354369", "AJ871581",
                       "KQ949024", "KP830094", "KF647220", "JXDG01000003",
                       "CP000667", "FN689524", "KB913022", "LN879412",
                       "CM000176", "KU946987", "HQ731031"])

        is_known = False
        if gbid in pos_set:
            is_known = True

        homologs_set = set(homologs_list)
        any_homologs_in_pos = bool(homologs_set & pos_set)

        if homologs == 0:
            homologs_list = None

        # Print only hilghly ranked clusters from the 616 novel targets list
        # if score_final < 3:
        #     continue

        f.write("%s\t%s\t%s\t%s\t" %
                (target_name, target_id, gbid, score_final))

        f.write("%s\t%s\t%s\t%s\t%s\t" % (d, homologs, homologs_list, target_count, copynum))
        f.write("%s\t%s\t%s\t" % (Sp_corr, p_val2, score))

        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
                (antismashfile, cluster_num, cluster_type, gene_name,
                 species, phyla, is_known, any_homologs_in_pos))

        f.write("%s\t%s\t%s\t%s\t%s\n" %
                (score_distance, score_homologs,
                 score_target_frequency, score_copy, score_coevolution))

    f.close()

if __name__ == "__main__":
    main()
