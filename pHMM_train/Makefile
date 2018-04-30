################################################################################################################################## 
# Compare 4 pHMMs interms of how they score:
# 1.  a positive training set (known TypeI, Type II, cisatpksnrps, transatpksnrps KS)
# 2.  a negative training set (FabF, TypeIII, TypeI, Type II, cisatpksnrps, transatpksnrps KS)
# Set of HMMs tested (for now):
# t1ks.11.hmm
# t2ks.12.hmm
# cisatpksnrps.10.hmm
# transatpksnrps.7.hmm
# Choose cutoff by testing HMMs against positive and negative train sets
# HMMs were made manually and stored in the dir HMM/

hmm_dir = HMM
# Positive training sets
fasta_dir = FastaTestSets

# Store results of hmmer search in these dirs
hmm_run_dir = $(hmm_dir).hmmsearch.$(fasta_dir)

# Store results of parsed hmmer results in these dirs.
# These parsed results can be read by R.
hmm_parse_dir = $(hmm_run_dir).parse

script_dir = scripts

# hmmsearch the selected hmms  vs. selected fasta files
$(hmm_run_dir): $(hmm_dir) $(fasta_dir)
	perl $(script_dir)/searchHMMs.pl $^ $@; \

$(hmm_parse_dir): $(hmm_run_dir) $(fasta_dir)
	perl $(script_dir)/parseHMMsScores.pl $^ $@; \

hmmsearch: $(hmm_run_dir) $(hmm_parse_dir)
