#!/usr/bin/env python
import sys
import tqdm
import subprocess
import os


gbids = set()

antismash_dir = "antismash_output_assemblies_all/"
antismash_dir2 = "antismash_output_assemblies_all_deleted/"
for gbidfull in tqdm.tqdm(os.listdir(antismash_dir)):
    

#    if not gbidfull.startswith("BDBD01000011"):
#        continue
       
#    print gbidfull

    if len(gbidfull.split("_")) < 2:
        gbids.add(gbidfull)


for gbidfull in tqdm.tqdm(os.listdir(antismash_dir)):
    
#    if not gbidfull.startswith("BDBD01000011"):
#        continue
       
#    print gbidfull

    if len(gbidfull.split("_")) >= 2:
        gbid = gbidfull.rsplit("_", 2)[0]
        blast_coord_start, blast_coord_end = gbidfull.split("_")[1:]
        if gbid in gbids:
            print 'DELETE %s gbid %s' % (gbid, gbidfull)
            subprocess.call(["mv", os.path.join(antismash_dir, gbid), os.path.join(antismash_dir2, gbid)]) 
            print ["mv", os.path.join(antismash_dir, gbid), os.path.join(antismash_dir2, gbid)]
            gbids.discard(gbid)

