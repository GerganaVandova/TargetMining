#!/usr/bin/env python
from Bio import SeqIO

f = open("targets.616.fa.cleannames", "w")

for record in SeqIO.parse(open("../Antismash_gbids/targets.616.fa.cleannames.withec", "rU"), "fasta"):
    # >DEG10180034_Glutamate-1-semialdehyde_2,1-aminomutase_(EC_5.4.3.8)
    gbidfull = record.id
    new_gbid = gbidfull.split("_(EC_")[0]
    seq = str(record.seq)
    f.write(">%s\n" % new_gbid)
    f.write(seq)
    f.write("\n")

f.close()
