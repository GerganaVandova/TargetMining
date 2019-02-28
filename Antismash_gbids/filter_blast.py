#!/usr/bin/env python


def main():

    E_VALUE_THRESH = float(1e-8)
    IDENTITY_THRESH_INIT = float(0.3)
    IDENTITY_THRESH_FABS = float(0.6)

    blastp_file = "out.92"
    f = open(blastp_file).readlines()

    #params = [blastp_file, str(E_VALUE_THRESH), str(IDENTITY_THRESH_INIT)]
    #blastp_filtered = ".".join(params)
    blastp_filtered = blastp_file + ".filtered"
    outf = open(blastp_filtered, "w")
    # qseqid sseqid sstart send nident qlen slen evalue
    for line in f:
        line = line.strip()
        print line
        features = line.split("\t")
        feats = str('\t'.join(features))
        qseqid, sseqid, sstart, send, nident, qlen, slen, evalue = features
        identity = float(nident)/float(qlen)

        #if "PtmP3_FabB-F" == qseqid:
        if "FAB" in qseqid:
            IDENTITY_THRESH = IDENTITY_THRESH_FABS
        else:
            IDENTITY_THRESH = IDENTITY_THRESH_INIT

        if float(evalue) < E_VALUE_THRESH and identity > IDENTITY_THRESH:
            print "%s\t%s\t%.2f\t%s" % (qseqid, sseqid, identity, evalue)
            outf.write("%s\t%s\t%.2f\t%s\n" %
                       (qseqid, sseqid, identity, evalue))
    outf.close()

if __name__ == "__main__":
    main()
