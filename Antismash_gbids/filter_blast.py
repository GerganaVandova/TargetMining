#!/usr/bin/env python


def main():

    E_VALUE_THRESH = float(1e-8)
    IDENTITY_THRESH_INIT = float(0.3)
    IDENTITY_THRESH_FABS = float(0.6)

    blastp_file = "out.616"
    f = open(blastp_file).readlines()

    blastp_filtered = blastp_file + ".filtered"
    outf = open(blastp_filtered, "w")
    # qseqid sseqid sstart send nident qlen slen evalue
    for line in f:
        line = line.strip()
        # print line
        features = line.split("\t")
        feats = str('\t'.join(features))
        qseqid, sseqid, sstart, send, nident, qlen, slen, evalue = features
        identity = float(nident)/float(qlen)

        if "FAB" in qseqid.upper():
            IDENTITY_THRESH = IDENTITY_THRESH_FABS
        else:
            IDENTITY_THRESH = IDENTITY_THRESH_INIT

        # For E. coli essential genes targets.616.fa
        fabs = ["DEG10180225", "DEG10180178", "DEG10180363", "DEG10180180", "DEG10180179"]
        if qseqid in fabs:
            IDENTITY_THRESH = IDENTITY_THRESH_FABS
        else:
            IDENTITY_THRESH = IDENTITY_THRESH_INIT
        
        #if "DEG10180225" in qseqid or "DEG10180178" in qseqid or "DEG10180363" in qseqid or "DEG10180180" in qseqid or "DEG10180179" in qseqid:
        #    IDENTITY_THRESH = IDENTITY_THRESH_FABS
        #else:
        #    IDENTITY_THRESH = IDENTITY_THRESH_INIT

        if float(evalue) < E_VALUE_THRESH and identity > IDENTITY_THRESH:
            print "%s\t%s\t%.2f\t%s" % (qseqid, sseqid, identity, evalue)
            outf.write("%s\t%s\t%.2f\t%s\n" %
                   (qseqid, sseqid, identity, evalue))
    outf.close()

if __name__ == "__main__":
    main()
