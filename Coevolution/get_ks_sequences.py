#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import defaultdict


def get_KS_sequences(blast_file):
    # Get KS coodinates and sequences from Blast results

    fastagbids_to_coord = defaultdict(list)
    # Head of blast fasta file
    # CENS01067892.1__80_1441_marine
    # CP000510___1361206_1362423__0_1_4559598_
    for record in SeqIO.parse(open(blast_file, "rU"), "fasta"):
        gbidfull = record.id
        seq = record.seq
        gbid1, rest = gbidfull.split("__", 1)
        gbid = gbid1.split(".")[0]
        parts = filter(lambda x: x, rest.split('_'))
        coord = "-".join(parts[:2])
        fastagbids_to_coord[gbid].append((coord, seq))

    return(fastagbids_to_coord)


def get_target_sequences(antismash_sequences):
    # Get coordinates of target genes from all parsed antismash sequences

    antismashgbids_to_targetseq = {}

    # >ACXX02000001_cluster-1_14512-116691_transatpks-nrps_31720-32586_Cpap_\
    # 3701_acetyl-CoA carboxylase, carboxyl transferase, beta subunit
    for record in SeqIO.parse(open(antismash_sequences, "rU"), "fasta"):
        # if "JXDG01000003" not in record.id:
        #     continue
        parts = record.id.split("_")
        # some ids start with  NZ_: NZ_RFFG01000027
        if parts[0] == "NZ":
            cluster = "_".join(parts[:2])
            target_gene_coord = parts[5]
        else:
            cluster, _, _, _, target_gene_coord = record.id.split("_")[:5]
        target_seq = record.seq
        key = cluster + "_" + target_gene_coord
        antismashgbids_to_targetseq[key] = target_seq

    return(antismashgbids_to_targetseq)


def get_targets(antismash_outfilename):
    # Get list of target gene names

    targets = set()
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target = line.split("\t")[0]
        targets.add(target)

    return(targets)


def get_antismash_clusters(gene_name, antismash_outfilename):
    # Get coordinates of KS and a target gene from Antismash_output results

    antismashgbids_to_kscoord = defaultdict(list)
    antismashgbids_to_targetcoord = defaultdict(list)
    antismashgbids_to_clustercoord = defaultdict(list)

    # ['Target', 'Cluster', 'Clusternum', 'Clustercoord', 'Type', 'Gene', \
    # 'sstart', 'send', 'nident', 'querylen', 'slen', 'pident', 'evalue', \
    # 'KS start', 'KS end', 'Distance', 'Cluster len', 'PKS_KS', 'PKS_KR', \
    # 'PKS_DH', 'PKS_ER', 'KS_AT', 'KS_ACP', 'AMP-binding', \
    # 'Condensation_LCL', 'Condensation_DCL', 'Condensation_Starter', 'PCP']
    # AdmT_ACC  ACXX02000001    1   14512-116691    transatpks-nrps\
    # Cpap_3701_acetyl-CoA    31720   32586   118.0   304.0   288 0.39\
    #    1e-67       38012   39229   5426102179  13  9   7   0
    outfile = open(antismash_outfilename).readlines()[1:]
    for line in outfile:
        line = line.strip()
        target, cluster, _, clustercoord, _, target_gene, target_gene_start, \
            target_gene_end, _, _, _, _, _, \
            ks_start, ks_end = line.split("\t")[:15]
        if target != gene_name:
            continue

        target_gene_coord = target_gene_start + "-" + target_gene_end
        ks_coord = ks_start + "-" + ks_end
        antismashgbids_to_kscoord[cluster].append(ks_coord)
        antismashgbids_to_targetcoord[cluster].append(target_gene_coord)
        antismashgbids_to_clustercoord[cluster].append(clustercoord)

    return([antismashgbids_to_kscoord,
            antismashgbids_to_targetcoord,
            antismashgbids_to_clustercoord])


def get_fasta_files(gene_name,
                    antismashgbids_to_kscoord,
                    antismashgbids_to_targetcoord,
                    antismashgbids_to_clustercoord):
    # Match antismash coord with blast coord to extract KS DNA sequence

    ks_outfilename = gene_name + ".KS" + ".fasta"
    f = open(ks_outfilename, "w")
    for gbid in sorted(antismashgbids_to_kscoord.keys()):
        aks_coord = antismashgbids_to_kscoord[gbid][0]
        ks_sequences = fastagbids_to_coord[gbid]
        clustercoord = antismashgbids_to_clustercoord[gbid][0]
        # if len(antismashgbids_to_clustercoord[gbid]) > 1:
        #     print gbid, "KScoord: ", antismashgbids_to_kscoord[gbid]
        #     print gbid, "clustercoord: ", antismashgbids_to_clustercoord[gbid]
        for ks in ks_sequences:
            ks_coord, ks_seq = ks
            if aks_coord == ks_coord:
                print "%s.%s.%s.KS.%s" % (gene_name, gbid, clustercoord, ks_coord)
                f.write(">%s.%s.%s.KS.%s\n" % (gene_name, gbid, clustercoord, ks_coord))
                f.write(str(ks_seq))
                f.write("\n")
    f.close()

    # Match antismash gbid and target coord to extract target gene DNA sequence
    target_outfilename = gene_name + ".fasta"
    ff = open(target_outfilename, "w")
    for gbid in sorted(antismashgbids_to_targetcoord.keys()):
        target_gene_coord = antismashgbids_to_targetcoord[gbid][0]
        # if len(antismashgbids_to_clustercoord[gbid]) > 1:
        #     print gbid, "targetcoord: ", antismashgbids_to_targetcoord[gbid]
        #     print gbid, "clustercoord: ", antismashgbids_to_clustercoord[gbid]
        clustercoord = antismashgbids_to_clustercoord[gbid][0]
        key = gbid + "_" + target_gene_coord
        target_sequence = antismashgbids_to_targetseq[key]
        print gene_name, gbid, clustercoord, target_gene_coord
        ff.write(">%s.%s.%s.%s\n" % (gene_name, gbid, clustercoord, target_gene_coord))
        ff.write(str(target_sequence))
        ff.write("\n")
    ff.close()


blast_file = \
    "../Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
antismash_outfilename = \
    "../Antismash_gbids/out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa"
antismash_sequences = "../Antismash_gbids/sequences.110k.coord.redundant.fasta"

antismashgbids_to_targetseq = get_target_sequences(antismash_sequences)
fastagbids_to_coord = get_KS_sequences(blast_file)

for target in get_targets(antismash_outfilename):
    antismashgbids_to_kscoord, \
    antismashgbids_to_targetcoord, \
    antismashgbids_to_clustercoord = \
        get_antismash_clusters(target, antismash_outfilename)
    get_fasta_files(target,
                    antismashgbids_to_kscoord,
                    antismashgbids_to_targetcoord,
                    antismashgbids_to_clustercoord)
