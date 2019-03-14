# Plotting phylogenetic trees in R using the APE package
# modified from Maureen Hillenmeyer
# Jan-2-2019

# Load the tree file
library(ape)
dir <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
# filename <- "KS.92.10kb.fasta.withFabF.cdhit.90.mafft.FastTree"
filename <- "KS.609.5kb.fasta.filtered.withFabF.cdhit.90.mafft.FastTree"

# Choose root sequence set
rootset <- "FabF"
file <- paste(dir, filename, sep="")
MyTree <- read.tree(file)

#############################################################
# # Highlight specific sequences
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename3 <- "mibig_refset.10.mod"
description_file3 = paste(dir1, filename3, sep="")
descriptions3 <- read.table(description_file3, sep="\t", as.is=T, row.names=1)
desc3.df <- data.frame(descriptions3)
tip.label3.df <- data.frame(MyTree$tip.label, row.names=1)
desc3.reordered <- desc3.df[rownames(tip.label3.df),]# This is the key step that matches the tree tip names to the external description file

# filename2 <- "KS.92.10kb.fasta.descr.species"
filename2 <- "KS.609.5kb.fasta.filtered.descr"
description_file2 = paste(dir1, filename2, sep="")
descriptions2 <- read.table(description_file2, sep="\t", as.is=T, row.names=1)
desc2.df <- data.frame(descriptions2)
tip.label2.df <- data.frame(MyTree$tip.label, row.names=1)
desc2.reordered <- desc2.df[rownames(tip.label2.df),]# This is the key step that matches the tree tip names to the external description file

# filename1 <- "KS.92.10kb.fasta.target"
filename1 <- "KS.609.5kb.fasta.filtered.target"
description_file = paste(dir1, filename1, sep="")
descriptions <- read.table(description_file, sep="\t", as.is=T, row.names=1)
desc.df <- data.frame(descriptions)
tip.label.df <- data.frame(MyTree$tip.label, row.names=1)
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
# Read the 2f file of fastaID and description
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
# filename1 <- "KS.92.10kb.fasta.phyla" # to get colors by phyla for KSs
filename1 <- "KS.609.5kb.fasta.filtered.phyla" # to get colors by phyla for KSs
phylum_file = paste(dir1, filename1, sep="")
phyla <- read.table(phylum_file, sep="\t", as.is=T, row.names=1)
phyla.df <- data.frame(phyla)
tip.label.df <- data.frame(MyTree$tip.label, row.names=1)
phyla.reordered <- phyla.df[rownames(tip.label.df),]  # reorder phyla

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

# Other bacteria 92 targets
acido <- phyla.reordered=="Acidobacteria" # gram neg
elusim <- phyla.reordered=="Elusimicrobia" # gram neg
fibro <- phyla.reordered=="Fibrobacteres" # gram neg?
gemma <- phyla.reordered=="Gemmatimonadetes" # gram neg
niro <- phyla.reordered=="Nitrospinae/Tectomicrobia group"
stram <- phyla.reordered=="Stramenopiles" # ?
verru <- phyla.reordered=="Verrucomicrobia" # gram neg
fungi <-phyla.reordered=="Fungi"
met <-phyla.reordered=="Metazoa"

# Other bacteria 609 targets
amo <-phyla.reordered=="Amoebozoa"
canc <-phyla.reordered=="Candidatus Cloacimonetes"
canr <-phyla.reordered=="Candidatus Riflebacteria"
chla <-phyla.reordered=="Chlamydiae"
chry <-phyla.reordered=="Chrysiogenetes"
eury <-phyla.reordered=="Euryarchaeota"
hapto <-phyla.reordered=="Haptophyceae"
nitrot <-phyla.reordered=="Nitrospinae/Tectomicrobia group"
nitro <-phyla.reordered=="Nitrospirae"

# colors
myCols1 <- c(rep("lightgray",length(MyTree$tip.label)))

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
#New colors 92 targets
myCols1[acido]="#acc6f6"
myCols1[elusim]="#a0db8e"
myCols1[fibro]="#cc5c5c"
myCols1[gemma]="#ff7373"
myCols1[niro]="#904690"
myCols1[stram]="#ffa500"
myCols1[verru]="#ffd5d5"
myCols1[fungi]="black"
myCols1[met]="black"
# New colors 609 targets
myCols1[amo]="grey"
myCols1[canc]="darkgreen"
myCols1[canr]="seagreen"
myCols1[chla]="palegreen"
myCols1[chry]="olivedrab"
myCols1[eury]="mediumspringgreen"
myCols1[hapto]="gold"
myCols1[nitro]="sandybrown"
myCols1[nitrot]="sandybrown"

myBG1 <- myCols1

##########################################################
# Plot rectangular phylogram With query sequences highlighted

outgroup <- grep(rootset, MyTree$tip.label, perl=TRUE)
MyTree.rooted <- root(MyTree,outgroup,node = NULL)
MyTree.ladderized <- ladderize(MyTree.rooted)

# mywidth=4; myheight=6 #for small rooted tree
# mywidth=6; myheight=22 #for 92 targets rooted tree
mywidth=20; myheight=60 #for 609 targets rooted tree

outfile <- paste(dir, filename, ".phyla.labeled.png", sep="")
title <- "KS.609targets.5kb.1166.sequences"
pdf(file=outfile, width=mywidth, height=myheight)
plot(MyTree.ladderized, main=title, font=.1, type="phylogram", edge.color="gray",
     edge.width=.5, show.tip.label=F, open.angle=5) # for rooted tree
# plot(MyTree, font=1, type="unrooted", edge.color="gray", edge.width=.5, show.tip.label=F, open.angle=5) # for unrooted tree
# tiplabels(pch=21, cex=.5, col=myCols, bg=myBG)# pch=21 circles, colored by target
tiplabels(pch=21, cex=.5, col=myCols1, bg=myBG1)# colored by phyla

# Print labels
# tiplabels(MyTree$tip.label, cex=.2, frame="none", adj=0) # full label
# tiplabels(phyla.reordered, cex=.1, frame="none", adj=0) # phyla label
# tiplabels(phyla.reordered, cex=.5, frame="none", adj=0) # phyla  label 609

tiplabels(desc2.reordered, cex=0.2, frame="none", adj=0) # short descr label
# tiplabels(desc3.reordered, cex=0.3, frame="none", adj=0) # to highlight mibig ref sequences

# nodelabels(MyTree$node.label, frame="none", cex=.1)
# edgelabels(MyTree$edge.label, frame="none", cex=.2)

add.scale.bar(cex=.5, lwd=.5)
dev.off()
