#################################################################
# Cat all files into one

script_dir = scripts
parse_blast_dir = blast_results_parse
blast_seqs_dir = blast_results_seqs
blast_results_tab_ks_file = $(blast_seqs_dir)/blast_results.parsed.KS.9f

$(blast_results_tab_ks_file):
	mkdir -p $(blast_seqs_dir)
	cat $(parse_blast_dir)/KS* > $@

concat: $(blast_results_tab_ks_file)

#################################################################
# Uniqify the proteins
# eliminiate 100% identical proteins here at this step.
# Choose longest nucleotide synonym.
# Also prioritize the Mibig set of nucleotide IDs:  Ie if this protein exists 
# in multiple nucleotide records, one of which is Mibig's nucleotide, 
# then keep the Mibig nucleotide, not the other one.
# Mark synonyms in a separate file: filename.synonyms
# length column is 6 (starting from 0) and seq col is 9
# * NOTE this script also removes dashes and asterisks


$(blast_results_tab_ks_file).longest_synonym: $(blast_results_tab_ks_file)
	
	perl $(script_dir)/chooseLongestSeqSynonym.pl $^ 6 8 $(blast_results_tab_ks_file).synonyms mibig.gbids.v1.4.pks > $@

longest: $(blast_results_tab_ks_file).longest_synonym $(blast_results_tab_ks_file).longest_synonym


# Compare the positive set (clusters chosen to be syntesized from prev time I ran the pipeline) to the found set - found them all?

$(blast_results_tab_ks_file).longest_synonym.compareposset: $(blast_results_tab_ks_file).longest_synonym positive_set.ids
	perl $(script_dir)/compareErin.pl $^ > $@

compare: $(blast_results_tab_ks_file).longest_synonym.compareposset


#################################################################
# Make fasta file of blast hit sequences to use in HMM search (without using min cluster size)

blast_ks_fasta = $(blast_seqs_dir)/blast_results.KS.fasta

$(blast_ks_fasta): $(blast_results_tab_ks_file).longest_synonym
	cat $^ \
	| $(script_dir)/multistab2fasta.pl $^ > $@

desclen2 = 200
blast_ks_fasta_cleanname = $(blast_ks_fasta).cleanName

$(blast_ks_fasta_cleanname): $(blast_ks_fasta)
	perl $(script_dir)/cleanName.pl $^ $(desclen2) > $@

fasta: $(blast_ks_fasta_cleanname)

clean_fasta: 
	rm -f $(blast_clf_fasta_cleanname) $(blast_ks_fasta_cleanname)
	rm -f $(blast_clf_fasta) $(blast_ks_fasta)


########################################################################################
# GV 031615 Eliminate identical KS/CLFs/ACPs from this set

cdhit_cutoff = .99
cdhit_bin = ~maureenh/Cdhit/cd-hit-v4.6.1-2012-08-27/cd-hit
name_length = 200

$(blast_ks_fasta).cleanName.cdhit$(cdhit_cutoff): $(blast_ks_fasta).cleanName
	$(cdhit_bin) -i $^ -o $@ -c $(cdhit_cutoff) -d $(name_length) -M 10000

cdhit: $(blast_ks_fasta).cleanName.cdhit$(cdhit_cutoff) $(blast_ks_fasta).cleanName.cdhit$(cdhit_cutoff)

all: cdhit compare
                                                                                                       
