#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

E_VALUE_THRESH = float(1e-8)
# Remove identity threshold for the 619 essential E coli genes
IDENTITY_THRESH = float(0.2)
IDENTITY_THRESH_FABS = float(0.2)

input_file = sys.argv[1]  # out.target.619.all

f = open(input_file).readlines()
filename_filtered = sys.argv[2]  # out.targets.619.filtered

ff = open(filename_filtered, "w")
# >NZ_ANBB01000037_cluster-1_75523-118144_t2fas-fabH_75248-75955_D459_RS33635_hypothetical protein

# 1_75523-118144_t2fas-fabH_75248-75955_D459_RS33635_hypothetical protein
# 1_79330-123547_t2fas-UbiA_cyclase_99330-100598_ctg1_orf00116_- cyclase

for line in f:
    qseqid, sseqid, sstart, send, nident, qlen, slen, evalue = line.split("\t")
    gbid, descr = sseqid.split("_cluster-")

    clusternum = descr.split("_")[0]
    coord = descr.split("_")[1]
    if "UbiA_cyclase" in descr:
        descr = descr.replace("UbiA_cyclase", "UbiA-cyclase")
    if "head_to_tail" in descr:
	descr = descr.replace("head_to_tail", "head-to-tail")
    clustertype = descr.split("_")[2]
    prot_coord = descr.split("_")[3]
    print "\nqseqid: ", qseqid, "\nsseqid: ", sseqid, "\nprot_coord: ", prot_coord
    prot_start, prot_end = prot_coord.split("-")
    protname = ("_").join(descr.split("_")[4:])

    nident = float(nident)
    qlen = float(qlen)
    identity = float(nident)/qlen
    evalue = float(evalue)
    if "Fab" in qseqid:
        if evalue < E_VALUE_THRESH and identity > IDENTITY_THRESH_FABS:
            # write ture protein coordinates
            print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s" % (qseqid,
                   gbid, clusternum, coord, clustertype, protname,
                   prot_start, prot_end, nident, qlen, slen, identity, evalue))
            ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s" % (qseqid,
                     gbid, clusternum, coord, clustertype, protname,
                     prot_start, prot_end, nident, qlen, slen, identity, evalue))
            ff.write("\n")
    elif evalue < E_VALUE_THRESH and identity > IDENTITY_THRESH:
            # identity threshold = 0.2 for the 619 essential E coli genes
            print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s" % (qseqid,
                   gbid, clusternum, coord, clustertype, protname,
                   prot_start, prot_end, nident, qlen, slen, identity, evalue))
            ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s" % (qseqid,
                     gbid, clusternum, coord, clustertype, protname,
                     prot_start, prot_end, nident, qlen, slen, identity, evalue))
            ff.write("\n")

ff.close()
