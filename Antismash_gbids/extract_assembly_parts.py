#!/usr/bin/env python
import subprocess
import os
import tqdm
import signal
import sys
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
from multiprocessing import Pool

TEMP_ANTISMASH_OUTPUT = "/home/gvandova/antismash_output_assemblies_2018/"
FINAL_ANTISMASH_OUTPUT = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Antismash_gbids/antismash_output_assemblies"
GB_FOLDER = "/home/gvandova/gbdir_assemblies"

LOOKOUT_RANGE = 150000

def select_ranges(ranges):
    final_ranges = []
    for x, y in ranges:
        x = max(0, x - LOOKOUT_RANGE)
        y = y + LOOKOUT_RANGE
        final_ranges.append((x, y))
    final_ranges = sorted(final_ranges)
    i = 0
    while i < len(final_ranges) - 1:
        x, y = final_ranges[i]
        xx, yy = final_ranges[i + 1]
        if y >= xx:
            min_x = min(x, xx)
            max_y = max(y, yy)
            final_ranges[i] = (min_x, max_y)
            final_ranges.pop(i + 1)
        else:
            i += 1

    return final_ranges

if __name__ == "__main__":
    # Get coordinates of KS from Blast results 
    fasta_file = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Blast/blast_results_seqs/blast_results.KS.fasta.cleanName.cdhit.99"
    gbids_to_coord = defaultdict(list)
    # CENS01067892.1__80_1441_marine
    # CP000510___1361206_1362423__0_1_4559598_
    for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
        gbidfull = record.id
        print gbidfull
        gbid, rest = gbidfull.split("__", 1)
        parts = filter(lambda x: x, rest.split('_'))
        start = int(parts[0])
        end = int(parts[1])
        print gbid, start, end
        gbids_to_coord[gbid].append((start, end))

    gbdir = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/assembly_gb/"
    gbids = set()
    gbidfile = "/mnt/gnpn/gnpn/projects/orphanpks/TargetMining/Genbank/gbids.assembly.unique.txt"
    ff = open(gbidfile, "r")
    for line in ff:
        gbid = line.strip()
        gbids.add(gbid)

    count = 0

    for gbid in tqdm.tqdm(sorted(gbids_to_coord.keys())):
        if gbid not in gbids:
            continue

        ranges = select_ranges(gbids_to_coord[gbid])
        min_coord = min([x for x, y in gbids_to_coord[gbid]])
        max_coord = max([y for x, y in gbids_to_coord[gbid]])

        #count += 1
        assembly_gbfile = os.path.join(gbdir, gbid + ".gbff")
        #print count, gbid, start, end, assembly_gbfile

        records = list(SeqIO.parse(open(assembly_gbfile, "rU"), "genbank"))
        if len(records) > 1:
            print "ERROR: expecting single seq %s" % assembly_gbfile
            sys.exit(-1)


        r = records[0]
        seqlen = len(r.seq)
        if seqlen < 1000000:
            continue
        
        ranges = select_ranges(gbids_to_coord[gbid])
        min_coord = min([x for x, y in gbids_to_coord[gbid]])
        max_coord = max([y for x, y in gbids_to_coord[gbid]])


#        for x, y in gbids_to_coord[gbid]:
#            print x, y
#        print gbid, min_coord, max_coord
#        print "Final", gbid, min_coord, max_coord, ranges

        count += 1
#        print "BIG %s %s %s" % (gbid, min_coord, max_coord)
        i = 0
        for x, y in ranges:
            y = min(y, seqlen)
            s = SeqRecord(r.seq[x:y], r.id, "", "")
            outfile = os.path.join(gbdir, "cut_seq", "%s_%s_%s.gbff" % (gbid, x, y))
            SeqIO.write(s, outfile, "genbank")
#            sys.exit(0)
            i += 1

        #print count, gbid, start, end, assembly_gbfile, seqlen
            #    print seq
            #SeqIO.write(record, output_file, 'genbank')
            #outf.write(gbid)
            #outf.write("\t")
            #outf.write(gbid_to_filename[gbid])
            #outf.write("\n")

    print count

    sys.exit(0)
    # Run Antismash 
    total = len(gbids)
    print "Total %s" % total
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    p = Pool(processes=45)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        for _ in tqdm.tqdm(p.imap_unordered(f, gbids), total=len(gbids)):
            pass

    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        p.terminate()
    else:
        print("Normal termination")
        p.close()
    p.join()

    print "Done, exiting"

