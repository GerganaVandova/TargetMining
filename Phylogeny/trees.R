# Plotting phylogenetic trees in R using the APE package
# modified from Maureen Hillenmeyer
# Jan-2-2019

# Load the tree file
library(ape)
dir <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename <- "KS.92.5kb.fasta.withFabF.cdhit.90.mafft.FastTree"

# Choose root sequence set
rootset <- "FabF"

file <- paste(dir, filename, sep="")
MyTree <- read.tree(file)

#############################################################
# # Highlight specific sequences
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
#filename1 <-"out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.cluster_type"
# whole description in second filed:
# ACXX02000001_38012-39229	ACXX02000001_38012-39229_AdmT_ACC_transatpks-nrps

# filename1 <- "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.descr"
filename3 <- "mibig_refset.10.mod"
description_file3 = paste(dir1, filename3, sep="")
descriptions3 <- read.table(description_file3, sep="\t", as.is=T, row.names=1)
desc3.df <- data.frame(descriptions3)
tip.label3.df <- data.frame(MyTree$tip.label, row.names=1)
# reorder descriptions
desc3.reordered <- desc3.df[rownames(tip.label3.df),]# This is the key step that matches the tree tip names to the external description file

# filename2 <- "KS.12.10kb.fasta.descr.species"
filename2 <- "KS.92.5kb.fasta.descr.species"
description_file2 = paste(dir1, filename2, sep="")
descriptions2 <- read.table(description_file2, sep="\t", as.is=T, row.names=1)
desc2.df <- data.frame(descriptions2)
tip.label2.df <- data.frame(MyTree$tip.label, row.names=1)
# reorder descriptions
desc2.reordered <- desc2.df[rownames(tip.label2.df),]# This is the key step that matches the tree tip names to the external description file
# Assign targets to variables


# filename1 <- "KS.12.10kb.fasta.target"
filename1 <- "KS.92.5kb.fasta.target"
description_file = paste(dir1, filename1, sep="")
descriptions <- read.table(description_file, sep="\t", as.is=T, row.names=1)
desc.df <- data.frame(descriptions)
tip.label.df <- data.frame(MyTree$tip.label, row.names=1)
# reorder descriptions
desc.reordered <- desc.df[rownames(tip.label.df),]# This is the key step that matches the tree tip names to the external description file
# Assign targets to variables
admt <- desc.reordered=="AdmT_ACC"
sal <- desc.reordered=="SalI_beta_proteasome"
dnan <- desc.reordered=="GriR_DnaN"
eftu <- desc.reordered=="EF-Tu"
fabb <- desc.reordered=="PtmP3_FabB-F"
fabi <- desc.reordered=="BatG_FabI"
gyrb <- desc.reordered=="GyrB-R"
ile <- desc.reordered=="mupM_Ile-tRNA-syn"
thr <- desc.reordered=="borI_Thr-tRNA-syn"
leu <- desc.reordered=="agnB2_Leu-tRNA-syn"
rub <- desc.reordered=="rubR1_TIF"
trp <- desc.reordered=="Ind0_Trp-tRNA-syn"
# colors
myCols <- c(rep("black",length(MyTree$tip.label)))
myCols[admt]="blue"
myCols[sal]="lightblue"
myCols[dnan]="cyan"
myCols[eftu]="darkblue"
myCols[fabb]="red"
myCols[fabi]="purple"
myCols[gyrb]="magenta"
myCols[ile]="brown"
myCols[thr]="lightpink"
myCols[leu]="orange"
myCols[rub]="green"
myCols[trp]="lightgreen"

myBG <- myCols

##############################################################
#Descriptions

# Read the 2f file of fastaID and description
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename1 <- "KS.92.5kb.fasta.phyla" # to get colors by phyla for KSs
phylum_file = paste(dir1, filename1, sep="")
phyla <- read.table(phylum_file, sep="\t", as.is=T, row.names=1)
phyla.df <- data.frame(phyla)
tip.label.df <- data.frame(MyTree$tip.label, row.names=1)
# reorder phyla
phyla.reordered <- phyla.df[rownames(tip.label.df),]

# Identify various bacterial phyla in the phylum file, assign to variables
actino <- phyla.reordered=="Actinobacteria"
proteo <- phyla.reordered=="Proteobacteria"
firm <- phyla.reordered=="Firmicutes"
bactero <- phyla.reordered=="Bacteroidetes"
cyano <- phyla.reordered=="Cyanobacteria"
spiro <- phyla.reordered=="Spirochaetes"
verru <- phyla.reordered=="Verrucomicrobia"
plancto <- phyla.reordered=="Planctomycetes"
thermo <- phyla.reordered=="Thermotogae"
chloroflexi <- phyla.reordered=="Chloroflexi"
syner <- phyla.reordered=="Synergistetes"
aqui <- phyla.reordered=="Aquificae"
cloa <- phyla.reordered=="Cloacimonetes"
# colors
myCols1 <- c(rep("black",length(MyTree$tip.label)))

# Terrabacteria: blues (per Hedges MBE 2009)
myCols1[actino]="blue"
myCols1[firm]="lightblue"## old: "orange"
myCols1[cyano]="cyan"
myCols1[chloroflexi]="darkblue"
# Hydrobacteria (per Hedges MBE 2009)
myCols1[proteo]="red"
myCols1[bactero]="purple"
myCols1[plancto]="magenta"
myCols1[verru]="brown"
myCols1[spiro]="lightpink"
# Other bacteria (per Hedges MBE 2009)
myCols1[thermo]="orange"
myCols1[aqui]="orange"
myCols1[syner]="orange"
myCols1[cloa]="orange"
# Other bacteria (random, check recommendations)

myBG1 <- myCols1

##########################################################
# Plot rectangular phylogram With query sequences highlighted

treetype <- "phylogram";
# treetype <- "unrooted";
outgroup <- grep(rootset, MyTree$tip.label, perl=TRUE)
MyTree.rooted <- root(MyTree,outgroup,node = NULL)
MyTree.ladderized <- ladderize(MyTree.rooted)

mywidth=4; myheight=6 #for small rooted tree
# mywidth=10; myheight=30 #for big rooted tree
# mywidth=6; myheight=4 #for unrooted tree

edge.color <- "gray"
myPch <- 21 # circles

outfile <- paste(dir, filename, ".", treetype,".phyla.png", sep="")
pdf(file=outfile, width=mywidth, height=myheight)
plot(MyTree.ladderized, font=1, type=treetype, edge.color=edge.color, edge.width=.5, show.tip.label=F, open.angle=5) # for rooted tree
# plot(MyTree, font=1, type=treetype, edge.color=edge.color, edge.width=.5, show.tip.label=F, open.angle=5) # for unrooted tree

# colored dots
# selectPchCex=.5 # for coloring dots by phyla, no labels
selectPchCex=.15 # for coloring dots by phyla with labels
# tiplabels(pch=myPch, cex=selectPchCex, col=myCols, bg=myBG)# colored by target
tiplabels(pch=myPch, cex=selectPchCex, col=myCols1, bg=myBG1)# colored by phyla

# Print labels
selectCex <- 1
# tiplabels(MyTree$tip.label[printlabels], printlabels, cex=selectCex, frame="none", adj=0) ### comment out if don't want to show labels
allLabCex <- .5
# allLabCex <- .1
# tiplabels(MyTree$tip.label, cex=.1, frame="none", adj=0) ### comment out if don't want to show labels
# tiplabels(MyTree$tip.label, cex=.2, frame="none", adj=0) ### comment out if don't want to show labels
# nodelabels(MyTree$node.label, frame="none", cex=.1)
#edgelabels(MyTree$edge.label, frame="none", cex=.2)

# tiplabels(phyla.reordered, cex=.1, frame="none", adj=0) # if you want phyla displayed
tiplabels(desc2.reordered, cex=0.1, frame="none", adj=0) # to label short descr
# tiplabels(desc3.reordered, cex=0.3, frame="none", adj=0) # to highlight mibig ref sequences

add.scale.bar(cex=allLabCex, lwd=selectPchCex)
dev.off()
