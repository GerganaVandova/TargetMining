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

a1<- desc.reordered=='sp_P17109_MEND_ECOLI'
a2<- desc.reordered=='sp_P0A6K3_DEF_ECOLI'
a3<- desc.reordered=='mfR6_squalene_synthase'
a4<- desc.reordered=='sp_P16659_SYP_ECOLI'
a5<- desc.reordered=='sp_P05041_PABB_ECOLI'
a6<- desc.reordered=='sp_P0A8M3_SYT_ECOLI'
a7<- desc.reordered=='sp_P0A6W3_MRAY_ECOLI'
a8<- desc.reordered=='sp_P0A7Z4_RPOA_ECOLI'
a9<- desc.reordered=='sp_P60785_LEPA_ECOLI'
a10<- desc.reordered=='sp_P10443_DPO3A_ECOLI'
a11<- desc.reordered=='sp_Q9RDT5_WALR_STAAU'
a12<- desc.reordered=='sp_P21889_SYD_ECOLI'
a13<- desc.reordered=='sp_P0A6G7_CLPP_ECOLI'
a14<- desc.reordered=='tr_E2QIY7_E2QIY7_ECOLX'
a15<- desc.reordered=='sp_P0A6M8_EFG_ECOLI'
a16<- desc.reordered=='sp_P77781_MENI_ECOLI'
a17<- desc.reordered=='sp_P0A8L1_SYS_ECOLI'
a18<- desc.reordered=='sp_P04079_GUAA_ECOLI'
a19<- desc.reordered=='sp_P07862_DDLB_ECOLI'
a20<- desc.reordered=='mpaF_IMDH'
a21<- desc.reordered=='sp_P0CE47_EFTU1_ECOLI'
a22<- desc.reordered=='sp_P0AGA2_SECY_ECOLI'
a23<- desc.reordered=='sp_P28305_PABC_ECOLI'
a24<- desc.reordered=='sp_P0A725_LPXC_ECOLI'
a25<- desc.reordered=='sp_P10408_SECA_ECOLI'
a26<- desc.reordered=='sp_P0AG30_RHO_ECOLI'
a27<- desc.reordered=='sp_P0AGB6_RPOE_ECOLI'
a28<- desc.reordered=='sp_P08312_SYFA_ECOLI'
a29<- desc.reordered=='sp_P0AEK2_FABG_ECOLI'
a30<- desc.reordered=='sp_Q9RDT3_WALK_STAAU'
a31<- desc.reordered=='sp_P0A749_MURA_ECOLI'
a32<- desc.reordered=='sp_P45568_DXR_ECOLI'
a33<- desc.reordered=='sp_P0A9M0_LON_ECOLI'
a34<- desc.reordered=='sp_P08373_MURB_ECOLI'
a35<- desc.reordered=='sp_P0A6B4_ALR1_ECOLI'
a36<- desc.reordered=='beta_lactamase'
a37<- desc.reordered=='sp_P17169_GLMS_ECOLI'
a38<- desc.reordered=='sp_P06986_HIS8_ECOLI'
a39<- desc.reordered=='sp_P0AGJ9_SYY_ECOLI'
a40<- desc.reordered=='sp_P0AAI3_FTSH_ECOLI'
a41<- desc.reordered=='sp_Q75R59_NQRF_VIBAN'
a42<- desc.reordered=='sp_P0A6N4_EFP_ECOLI'
a43<- desc.reordered=='sp_P07003_POXB_ECOLI'
a44<- desc.reordered=='sp_P00903_PABA_ECOLI'
a45<- desc.reordered=='sp_P32166_MENA_ECOLI'
a46<- desc.reordered=='sp_P0ABU0_MENB_ECOLI'
a47<- desc.reordered=='sp_P0C0V0_DEGP_ECOLI'
a48<- desc.reordered=='sp_P37353_MENE_ECOLI'
a49<- desc.reordered=='sp_P03007_DPO3E_ECOLI'
a50<- desc.reordered=='sp_P18335_ARGD_ECOLI'

myCols[a1]="#FFFF00"
myCols[a2]="#1CE6FF"
myCols[a3]="#FF34FF"
myCols[a4]="#FF4A46"
myCols[a5]="#008941"
myCols[a6]="#006FA6"
myCols[a7]="#A30059"
myCols[a8]="#FFDBE5"
myCols[a9]="#7A4900"
myCols[a10]="#0000A6"
myCols[a11]="#63FFAC"
myCols[a12]="#B79762"
myCols[a13]="#004D43"
myCols[a14]="#8FB0FF"
myCols[a15]="#997D87"
myCols[a16]="#5A0007"
myCols[a17]="#809693"
myCols[a18]="#FEFFE6"
myCols[a19]="#1B4400"
myCols[a20]="#4FC601"
myCols[a21]="#3B5DFF"
myCols[a22]="#4A3B53"
myCols[a23]="#FF2F80"
myCols[a24]="#61615A"
myCols[a25]="#BA0900"
myCols[a26]="#6B7900"
myCols[a27]="#00C2A0"
myCols[a28]="#FFAA92"
myCols[a29]="#FF90C9"
myCols[a30]="#B903AA"
myCols[a31]="#D16100"
myCols[a32]="#DDEFFF"
myCols[a33]="#000035"
myCols[a34]="#7B4F4B"
myCols[a35]="#A1C299"
myCols[a36]="#300018"
myCols[a37]="#0AA6D8"
myCols[a38]="#013349"
myCols[a39]="#00846F"
myCols[a40]="#372101"
myCols[a41]="#FFB500"
myCols[a42]="#C2FFED"
myCols[a43]="#A079BF"
myCols[a44]="#CC0744"
myCols[a45]="#C0B9B2"
myCols[a46]="#C2FF99"
myCols[a47]="#001E09"
myCols[a48]="#00489C"
myCols[a49]="#6F0062"
myCols[a50]="#0CBD66"

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

outfile <- paste(dir, filename, ".", treetype,".targets.png", sep="")
pdf(file=outfile, width=mywidth, height=myheight)
plot(MyTree.ladderized, font=1, type=treetype, edge.color=edge.color, edge.width=.5, show.tip.label=F, open.angle=5) # for rooted tree
# plot(MyTree, font=1, type=treetype, edge.color=edge.color, edge.width=.5, show.tip.label=F, open.angle=5) # for unrooted tree

# colored dots
# selectPchCex=.5 # for coloring dots by phyla, no labels
selectPchCex=.15 # for coloring dots by phyla with labels
tiplabels(pch=myPch, cex=selectPchCex, col=myCols, bg=myBG)# colored by target
# tiplabels(pch=myPch, cex=selectPchCex, col=myCols1, bg=myBG1)# colored by phyla

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
