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


def get_coevolution_score_and_homologs(coevolutionfilename):
    gbid_to_score = defaultdict()
    gbid_to_homologs = defaultdict()
    coevolutionfile = open(coevolutionfilename).readlines()
    for line in coevolutionfile:
        # feats = line.split("\t")
        # if len(feats) < 6:
        gbid, target_id, count, score, homologs, homologs_list, pident_list, \
            redundant, redundant_list = line.split("\t")
        mytuple = (target_id, count, score)
        pident_list = pident_list.strip()
        redundant_list = redundant_list.strip()
        value = "|".join(mytuple)
        gbid_to_score[gbid] = value
        mytuple2 = (homologs, homologs_list, pident_list, redundant, redundant_list)
        value2 = "|".join(mytuple2)
        gbid_to_homologs[gbid] = value2

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
    # Get target name from targets fasta file
    # >AM889285|cluster-1|2454838-2501896|t1pks|2500895-2502196|ctg1_139
    gbid_to_gene = defaultdict()
    for record in SeqIO.parse(open(genes_fasta, "rU"), "fasta"):
        gbidfull = record.id
        gbid, _, _, _, gene_coord, gene_name = gbidfull.split("|")  # for E. coli 609 targets
        gene_start, gene_end = gene_coord.split("-")
        gbid_to_gene[gbid, gene_start, gene_end] = gene_name
    return gbid_to_gene


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
    gbid_to_score, gbid_to_homologs = get_coevolution_score_and_homologs(coevolutionfilename)

    second_copy_filename = "../Second_copy/out.second_copy.609.10kb.filtered"
    gbid_to_copy = get_second_copy(second_copy_filename)

    genes_fasta = "../Antismash_gbids/cluster_genes.all.fasta"
    gbid_to_gene = get_gene_name(genes_fasta)

    f = open("Clusters.609.10kb.txt", "w")
    # f.write("target_name\tgbid\tscore_final\tantismashfile\tcluster_num\tcluster_type\thomologs\thomologs_list\tpident_list\tredundant\tredundant_list\tks_start\tks_end\tgene_name\ttarget_start\ttarget_end\tcluster_coord\tspecies\tphyla\td\tscore_distance\tcount\tscore_target_frequency\tscore\tscore_coevolution\tcopynum\tcomplete\tscore_copy\tscore_homologs\tpositive\tpos_homolog\tpos_redundant\n")
    f.write("target_name\ttarget id\tgbid\tscore_final\tantismashfile\tcluster_num\tcluster_type\tdistance\tspecies\tphyla\thomologs\thomologs_list\tredundant\tredundant_list\tcount\tcopynum\tcoevolution score\tpositive\tpos_homolog\tpos_redundant\n")

    # Read KS fasta file and write phyla for each gbid

    for record in tqdm.tqdm(SeqIO.parse(open("KS.609.10kb.fasta", "rU"), "fasta")):
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
        # if target_id != "PtmP3_FabB-F":
        #     continue

        # if gbid != "CP025542":
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

        if gbid not in gbid_to_score.keys():
            coevolution_score = target_id + '|0|' + '-1'
        else:
            coevolution_score = gbid_to_score[gbid]
        target_ids, count, score = coevolution_score.split("|")
        score = float(score)
        count = int(count)

        if gbid not in gbid_to_homologs.keys():
            homologs = -1
            homologs_list = []
            pident_list = []
            redundant = -1
            redundant_list = []
        else:
            homologs, homologs_list, pident_list, redundant, redundant_list = \
                gbid_to_homologs[gbid].split("|")
            homologs = int(homologs)
            homologs_list = json.loads(homologs_list)
            redundant = int(redundant)
            redundant_list = json.loads(redundant_list)

        second_copy = gbid_to_copy[gbid]
        copynum, genome_size, complete = second_copy.split("|")
        copynum = int(copynum)

        # 1. Calclulate KS-target distance score
        score_distance = -10000
        # print "distance is ", d, type(d)
        if d <= 5000:
            score_distance = 20
        elif d > 5000 and d <= 10000:
            score_distance = 10
        elif d > 10000 and d <= 20000:
            score_distance = 5
        else:
            score_distance = 0
        # print score_distance

        # 2. Calculate score based on occurences of targets
        score_target_frequency = -1000
        if count > 1 and count <= 20:
            score_target_frequency = 10
        elif count > 20:
            score_target_frequency = 0
        else:
            score_target_frequency = 0

        # 3. Calclulate score based on KS-target coevolution
        score_coevolution = -1000
        if score > 0 and score <= 10:
            score_coevolution = 20
        elif score > 10 and score <= 20:
            score_coevolution = 10
        elif score > 20:
            score_coevolution = 5
        else:
            score_coevolution = 0

        # 4. Calclulate score based on second copy of target
        score_copy = -1000
        if copynum > 1:
            score_copy = 10
        elif copynum == 1 and complete == "False":
            score_copy = 5
        else:
            score_copy = 0

        # 5. Calclulate score based on number of homologs >80% but <90% identity
        score_homologs = -1000
        if homologs >= 1:
            score_homologs = 10
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

        # print homologs_list, redundant_list
        homologs_set = set(homologs_list)
        redundant_set = set(redundant_list)
        # print pos_set, homologs_set, redundant_set
        any_homologs_in_pos = bool(homologs_set & pos_set)
        any_redundant_in_pos = bool(redundant_set & pos_set)
        # print any_homologs_in_pos, any_redundant_in_pos

        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
                (target_name, target_id, gbid, score_final,
                 antismashfile, cluster_num, cluster_type, d,
                 species, phyla, homologs, homologs_list,
                 redundant, redundant_list, count, copynum, score))

        # f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
        #         (target_name, target_id, gbid, score_final,
        #          antismashfile, cluster_num, cluster_type,
        #          homologs, homologs_list, pident_list,
        #          redundant, redundant_list))
        #
        # f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
        #         (ks_start, ks_end,
        #          target_start, target_end, gene_name,
        #          cluster_coord,
        #          species, phyla))
        #
        # f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
        #         (d, score_distance,
        #          count, score_target_frequency,
        #          score, score_coevolution,
        #          copynum, complete, score_copy,
        #          score_homologs))
        #
        f.write("%s\t%s\t%s\n" % (is_known, any_homologs_in_pos, any_redundant_in_pos))

    f.close()

if __name__ == "__main__":
    main()
