# Plotting phylogenetic trees in R using the APE package
# Maureen Hillenmeyer
# Jan-2-2019

# Load the tree file
library(ape)
dir <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
#filename <- "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.withFabF.fasta.mafft.FastTree"

filename <- "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.withFabF.fasta.cdhit.90.mafft.FastTree"

# Choose root sequence set
rootset <- "FabF"

file <- paste(dir, filename, sep="")
MyTree <- read.tree(file)

#############################################################
# # Highlight specific sequences
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
#filename1 <- "mibig_refset.8" # to get Erin's ref set cluster names for KSs
#filename1 <-"out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.cluster_type"

# whole description in second filed:
# ACXX02000001_38012-39229	ACXX02000001_38012-39229_AdmT_ACC_transatpks-nrps
filename1 <- "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.descr"

description_file = paste(dir1, filename1, sep="")
descriptions <- read.table(description_file, sep="\t", as.is=T, row.names=1)
desc.df <- data.frame(descriptions)
tip.label.df <- data.frame(MyTree$tip.label, row.names=1)
# reorder descriptions
desc.reordered <- desc.df[rownames(tip.label.df),]# This is the key step that matches the tree tip names to the external description file
##############################################################
#Descriptions

# Read the 2f file of fastaID and description
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename1 <- "out.targets.12.eval.1e-8.pident.30.filtered.10000.allpks.domains.268.taxa.KS.fasta.phyla" # to get colors by phyla for KSs

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
myCols <- c(rep("black",length(MyTree$tip.label)))

# Terrabacteria: blues (per Hedges MBE 2009)
myCols[actino]="blue"
myCols[firm]="lightblue"## old: "orange"
myCols[cyano]="cyan"
myCols[chloroflexi]="darkblue"
# Hydrobacteria (per Hedges MBE 2009)
myCols[proteo]="red"
myCols[bactero]="purple"
myCols[plancto]="magenta"
myCols[verru]="brown"
myCols[spiro]="lightpink"
# Other bacteria (per Hedges MBE 2009)
myCols[thermo]="orange"
myCols[aqui]="orange"
myCols[syner]="orange"
myCols[cloa]="orange"
# Other bacteria (random, check recommendations)

myBG <- myCols

##########################################################
# colors
# for now just b&w
#myCols <- c(rep("gray",length(MyTree$tip.label)))
#myBG <- myCols
#myBG<- c(rep("white",length(MyTree$tip.label)))

##########################################################
# Plot rectangular phylogram
# With query sequences highlighted
# Must choose outgroup sequence to root on

treetype <- "phylogram";
outgroup <- grep(rootset, MyTree$tip.label, perl=TRUE)

MyTree.rooted <- root(MyTree,outgroup,node = NULL)
MyTree.ladderized <- ladderize(MyTree.rooted)

mywidth=3; myheight=6 #for small tree
# mywidth=10; myheight=30 #for big tree
#mywidth=40; myheight=100# really big for browsing names

edge.color <- "gray"
myPch <- 21# circles

outfile <- paste(dir, filename, ".",treetype,".png", sep="")
pdf(file=outfile, width=mywidth, height=myheight)
plot(MyTree.ladderized, font=1, type=treetype, edge.color=edge.color, edge.width=.5, show.tip.label=F, open.angle=5)

# colored dots
# selectPchCex=.5 # for big pdf (10.30)
selectPchCex=.15 #for small pdf (2.6)
tiplabels(pch=myPch, cex=selectPchCex, col=myCols, bg=myBG)# colored label dots

# Print select labels
#  identified using grep above (Erin set)
selectCex <- 1
# tiplabels(MyTree$tip.label[printlabels], printlabels, cex=selectCex, frame="none", adj=0) ### comment out if don't want to show labels

# Print all labels
allLabCex <- .5
# allLabCex <- .1
# tiplabels(MyTree$tip.label, cex=.1, frame="none", adj=0) ### comment out if don't want to show labels
# tiplabels(MyTree$tip.label, cex=.2, frame="none", adj=0) ### comment out if don't want to show labels
# nodelabels(MyTree$node.label, frame="none", cex=.1)
#edgelabels(MyTree$edge.label, frame="none", cex=.2)

tiplabels(desc.reordered, cex=0.1, frame="none", adj=0) # to highlight ref sequences when pdf(3,6)
# # tiplabels(desc.reordered, cex=.3, frame="none", adj=0) # to highlight ref sequences for big pdf (10.30)
# # tiplabels(phyla.reordered, cex=.1, frame="none", adj=0) # if you want phyla displayed

add.scale.bar(cex=allLabCex, lwd=selectPchCex)
dev.off()
