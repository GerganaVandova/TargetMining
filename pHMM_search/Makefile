#################################################################
# HMM search for new t2clfs
#
# Note that HMMs were made using hmmbuild manually in HMM/
# See HMM/README
#

# input blast files
blast_seqs_dir = ../Blast/blast_results_seqs
blast_ks_fasta = $(blast_seqs_dir)/blast_results.KS.fasta

blast_ks_fasta_cleanname = $(blast_ks_fasta).cleanName

# KS
hmm_model_filename_ks = t1ks.10.hmm
hmm_model_file_ks = ../pHMM_train/HMM/$(hmm_model_filename_ks)
hmm_out_file_ks = $(hmm_model_filename_ks).hmmsearch



# hmmsearch score cutoffs
hmm_score_cutoff_ks = 200 #based on ? KS type1

script_dir = scripts

# Run hmmsearch
$(hmm_out_file_ks): $(blast_ks_fasta_cleanname)
	hmmsearch $(hmm_model_file_ks) $^ > $@

run_hmmsearch: $(hmm_out_file_ks)

# Parse output
hmm_parsed_file_ks = $(hmm_model_filename_ks).hmmsearch.parsed

$(hmm_parsed_file_ks): $(hmm_out_file_ks)
	perl $(script_dir)/parseHMM.pl $^ > $@

hmmsearch_parse: $(hmm_parsed_file_ks)

# Filter results according to selected hmmsearch cutoff score
hmm_filtered_file_ks = $(hmm_model_filename_ks).hmmsearch.parsed.$(hmm_score_cutoff_ks)

$(hmm_filtered_file_ks): $(hmm_parsed_file_ks)
	cat $^ \
	| perl -anF/"\t"/ -e 'print "$$_" if ($$F[8] > $(hmm_score_cutoff_ks));'\
	> $@; \


filter: $(hmm_filtered_file_ks)

clean:
	rm $(hmm_out_file_ks)

#################################################################
# Fasta seqs from HMM results using to count unique sequences
hmm_filtered_file_ks_fasta = $(hmm_filtered_file_ks).fasta

# KS fasta
$(hmm_filtered_file_ks_fasta): $(hmm_filtered_file_ks)
	cat $^ \
	| $(script_dir)/multistab2fasta.pl $^ > $@


#add this line to see if it is going to remove the circular dependnecy
fasta: $(hmm_filtered_file_ks_fasta)

#CleanSeq (remove dashes introduced during the hmm searcg)
name_length = 200

#$(hmm_filtered_file_ks).cleanSeq.fasta: t1ks.10.hmm.hmmsearch.parsed.200.fasta
$(hmm_filtered_file_ks).cleanSeq.fasta: $(hmm_filtered_file_ks).fasta
	perl $(script_dir)/cleanSeq.pl $^ > $@

cleanseq: $(hmm_filtered_file_ks).cleanSeq.fasta



# Use cdhit to eliminate identical sequences
cdhit_cutoff = .99
cdhit_bin = ~maureenh/Cdhit/cd-hit-v4.6.1-2012-08-27/cd-hit

$(hmm_filtered_file_ks).cdhit$(cdhit_cutoff): $(hmm_filtered_file_ks).cleanSeq.fasta
	$(cdhit_bin) -i $^ -o $@ -c $(cdhit_cutoff) -d $(name_length)

cdhit: $(hmm_filtered_file_ks).cdhit$(cdhit_cutoff)

all: run_hmmsearch hmmsearch_parse filter cdhit
