#!/usr/bin/python
import sys
import json


def main():

    # Append taxa to out file. To execute script:

    gbid_to_taxa = {}
    # BDBI01000023	Bacteria	Actinobacteria	Corynebacteriales	Nocardiaceae	Nocardia
    taxa_filename = "../Genbank/taxa.txt"
    taxa_file = open(taxa_filename).readlines()
    for line in taxa_file:
        line = line.strip()
        if len(line.split("\t", 1)) < 2:
            continue
        gbid, rest = line.split("\t", 1)
        taxa = "|".join(rest.split("\t"))
        gbid_to_taxa[gbid] = taxa

    input_file = open("out.12.filtered.10kb").readlines()
    output_file = "out.12.filtered.10kb.taxa"
    outf = open(output_file, "w")

    # ACXX02000001|AdmT_ACC|37972|38377|31720|32586|cluster-1|transatpks-nrps|14512-116691|5386
    for line in input_file:
        line = line.strip()
        gbid, rest = line.split("|", 1)
        if gbid not in gbid_to_taxa.keys():
            outf.write("%s|None\n" % line)
            print "No taxa", gbid
            continue
        print "%s|%s" % (line, gbid_to_taxa[gbid])
        outf.write("%s|%s\n" % (line, gbid_to_taxa[gbid]))

    outf.close()

if __name__ == "__main__":
    main()
