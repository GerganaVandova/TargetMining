# GV April-30-2018 Function to plot scores, given one pHMM,
# Assumes that pdf() has been called prior to function call

# Run script:
# /Library/Frameworks/R.framework/Versions/3.1/Resources/bin/Rscript selectCutoff.R


plotScores <- function(hmm_filename, name_text, mylabs, mylabs_readable,
						cutoff, yaxis) {

	# Test sets
	dir <- "HMM.hmmsearch.FastaTestSets.parse/"

	x.increment <- 1/length(mylabs)
	xloc <- seq(x.increment,1, by=x.increment)
	# name <- paste(name_text, hmm_filename, sep="") # Name of plot has the name of hmm as well
	name <- paste(name_text)

	for (i in 1:length(mylabs)) {
		fn <- paste(hmm_dir, dir, hmm_filename, ".", mylabs[i], ".fasta", sep="")
		t <- read.table(fn)
		x.i <- rep(xloc[i], length(t$V1))
		if ( i == 1) {
			plot(x.i, t$V1, ylim=c(0, yaxis), xlim=c(xloc[1], xloc[length(xloc)]), xaxt='n', xlab="", ylab="HMMER3 pHMM scores", main=name)
		} else {
			x.i <- rep(xloc[i], length(t$V1)); points(x.i, t$V1)
		}
	}

	# Axis labels
	mycex.axis=0.8 #before 0.8
	axis(side=1, at=xloc, labels=mylabs_readable, cex.axis=mycex.axis, lwd=0)#, padj=c(0,.5))
	abline(h=cutoff,col=2,lty=2)

}

hmm_dir <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/pHMM_train/"

# Print the output figure to pdf
outfile <- paste(hmm_dir,"compare.sets.0507.Fig2b.pdf", sep="")
pdf(file=outfile, width=8, height=12) #for plotting all panels
par(mfrow=c(4, 2)) # to plot all but KS and CLF domains

all_labs <- list(
	# c("TypeI.KS.5", "Cisat.KS.6", "Transat.KS.5", "TypeII.KS.30", "KSIII.10", "TypeIII.KS.3", "FabF"),
	c("TypeI.KS.5", "Cisat.KS.6", "Transat.KS.5", "TypeII.KS.30", "KSIII.10", "TypeIII.KS.3", "FabF"),
	c("Cisat.KS.6", "TypeI.KS.5", "Transat.KS.5", "TypeII.KS.30", "KSIII.10", "TypeIII.KS.3", "FabF"),
	c("Transat.KS.5", "TypeI.KS.5", "Cisat.KS.6", "TypeII.KS.30", "KSIII.10", "TypeIII.KS.3", "FabF"),
	c("TypeII.KS.30", "TypeI.KS.5", "Cisat.KS.6", "Transat.KS.5", "KSIII.10", "TypeIII.KS.3", "FabF")
)

all_labs_readable <- list(
	# c("TypeI KS", "Cisat KS", "Transat KS", "TypeII KS", "KS III", "TypeIII KS", "FabF"),
	c("TypeI KS", "Cisat KS", "Transat KS", "TypeII KS", "KS III", "TypeIII KS", "FabF"),
	c("Cisat KS", "TypeI KS", "Transat KS", "TypeII KS", "KS III", "TypeIII KS", "FabF"),
	c("Transat KS", "TypeI KS", "Cisat KS", "TypeII KS", "KS III", "TypeIII KS", "FabF"),
	c("TypeII KS", "TypeI KS", "Cisat KS", "Transat KS", "KS III", "TypeIII KS", "FabF")
)

cutoffs <- c(200, 200, 200, 200) # for all domains
# cutoffs <- c(200, 200, 200, 200, 200) # for all domains
# cutoffs <- c(200, 200) # for typeI and typeII pks

names <- c(
	"typeI KS pHMM (19 protein sequences)\n",
	# "typeI KS pHMM (10 protein sequences)\n",
	"cisatpksnrps KS pHMM (10 protein sequences)\n",
	"transatpksnrps KS pHMM (7 protein sequences)\n",
	"type II KS pHMM (12 protein sequences)\n"
)

hmm_filenames <- c(
	# "t1ks.19.hmm",
	"t1ks.10.hmm",
	"cisatpksnrps.10.hmm",
	"transatpksnrps.7.hmm",
	"t2ks.12.hmm"
)

yaxes <- c(800, 800, 800, 800) # for all domains
# yaxes <- c(800, 800, 800, 800, 800) # for all domains
# yaxes <- c(800, 800) # for all domains


for (i in 1:length(names)) {
	plotScores(hmm_filenames[i], names[i], all_labs[[i]], all_labs_readable[[i]], cutoffs[i], yaxes[i])

}

dev.off()
