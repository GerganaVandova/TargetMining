# Plotting phylogenetic trees in R using the APE package
# modified from Maureen Hillenmeyer
# Jan-2-2019

# Load the tree file
library(ape)
dir <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename <- "KS.14.10kb.fasta.cdhit.90t.withFabF.mafft.FastTree"


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

filename2 <- "KS.14.10kb.fasta.descr"
# filename2 <- "KS.616.10kb.fasta.filtered.descr"
# filename2 <- "KS.616.10kb.fasta.filtered.descr.noname"

description_file2 = paste(dir1, filename2, sep="")
descriptions2 <- read.table(description_file2, sep="\t", as.is=T, row.names=1)
desc2.df <- data.frame(descriptions2)
tip.label2.df <- data.frame(MyTree$tip.label, row.names=1)
desc2.reordered <- desc2.df[rownames(tip.label2.df),]# This is the key step that matches the tree tip names to the external description file

filename1 <- "KS.14.10kb.fasta.filtered.target"
# filename1 <- "KS.616.10kb.fasta.filtered.target"
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
trp <- desc.reordered=="Ind0_Trp-tRNA-syn"
metap <- desc.reordered=="methionine_aminopeptidase"
otcase <- desc.reordered=="ArgK_OTCase"
ser <- desc.reordered=="seryl-tRNA_synthetase"

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
myCols[trp]="lightgreen"
myCols[metap]="green"
myCols[otcase]="deeppink"
myCols[ser]="slateblue"

#
# a1<-desc.reordered=='sp_P00934_THRC_ECOLI'
# a2<-desc.reordered=='sp_P65636_ODP2_STAAN'
# a3<-desc.reordered=='sp_P60785_LEPA_ECOLI'
# a4<-desc.reordered=='DEG10180273'
# a5<-desc.reordered=='sp_P07623_METAS_ECOLI'
# a6<-desc.reordered=='DEG10180503'
# a7<-desc.reordered=='DEG10180500'
# a8<-desc.reordered=='DEG10180501'
# a9<-desc.reordered=='DEG10180504'
# a10<-desc.reordered=='DEG10180505'
# a11<-desc.reordered=='sp_P06986_HIS8_ECOLI'
# a12<-desc.reordered=='sp_P17443_MURG_ECOLI'
# a13<-desc.reordered=='DEG10180481'
# a14<-desc.reordered=='DEG10180480'
# a15<-desc.reordered=='sp_P0A749_MURA_ECOLI'
# a16<-desc.reordered=='DEG10180366'
# a17<-desc.reordered=='sp_P46022_MTGA_ECOLI'
# a18<-desc.reordered=='DEG10180470'
# a19<-desc.reordered=='DEG10180570'
# a20<-desc.reordered=='DEG10180004'
# a21<-desc.reordered=='sp_P00903_PABA_ECOLI'
# a22<-desc.reordered=='sp_P0A725_LPXC_ECOLI'
# a23<-desc.reordered=='DEG10180498'
# a24<-desc.reordered=='sp_P0A6K3_DEF_ECOLI'
# a25<-desc.reordered=='DEG10180491'
# a26<-desc.reordered=='DEG10180494'
# a27<-desc.reordered=='DEG10180495'
# a28<-desc.reordered=='sp_P05041_PABB_ECOLI'
# a29<-desc.reordered=='sp_P0AG99_SECG_ECOLI'
# a30<-desc.reordered=='DEG10180099'
# a31<-desc.reordered=='DEG10180318'
# a32<-desc.reordered=='sp_P07395_SYFB_ECOLI'
# a33<-desc.reordered=='DEG10180315'
# a34<-desc.reordered=='sp_P0A6G7_CLPP_ECOLI'
# a35<-desc.reordered=='sp_P45568_DXR_ECOLI'
# a36<-desc.reordered=='DEG10180569'
# a37<-desc.reordered=='rubR1_TIF'
# a38<-desc.reordered=='DEG10180406'
# a39<-desc.reordered=='sp_P26459_APPC_ECOLI'
# a40<-desc.reordered=='DEG10180155'
# a41<-desc.reordered=='sp_P77790_DDPX_ECOLI'
# a42<-desc.reordered=='sp_P10408_SECA_ECOLI'
# a43<-desc.reordered=='sp_P0AGB6_RPOE_ECOLI'
# a44<-desc.reordered=='sp_P08312_SYFA_ECOLI'
# a45<-desc.reordered=='sp_Q9RDT3_WALK_STAAU'
# a46<-desc.reordered=='sp_P0A7Z4_RPOA_ECOLI'
# a47<-desc.reordered=='sp_P0ABJ9_CYDA_ECOLI'
# a48<-desc.reordered=='sp_P17169_GLMS_ECOLI'
# a49<-desc.reordered=='sp_P12995_BIOA_ECOLI'
# a50<-desc.reordered=='sp_P0A6N4_EFP_ECOLI'
# a51<-desc.reordered=='DEG10180381'
# a52<-desc.reordered=='DEG10180418'
# a53<-desc.reordered=='sp_P37353_MENE_ECOLI'
# a54<-desc.reordered=='sp_P03007_DPO3E_ECOLI'
# a55<-desc.reordered=='sp_P18335_ARGD_ECOLI'
# a56<-desc.reordered=='sp_P17109_MEND_ECOLI'
# a57<-desc.reordered=='sp_P99084_DLDH_STAAN'
# a58<-desc.reordered=='sp_P10443_DPO3A_ECOLI'
# a59<-desc.reordered=='sp_P0A6M8_EFG_ECOLI'
# a60<-desc.reordered=='sp_P0A6H1_CLPX_ECOLI'
# a61<-desc.reordered=='DEG10180039'
# a62<-desc.reordered=='sp_P0A6B4_ALR1_ECOLI'
# a63<-desc.reordered=='sp_P07862_DDLB_ECOLI'
# a64<-desc.reordered=='sp_P0CE47_EFTU1_ECOLI'
# a65<-desc.reordered=='sp_P0A6J8_DDLA_ECOLI'
# a66<-desc.reordered=='sp_P08373_MURB_ECOLI'
# a67<-desc.reordered=='sp_P9WN39_GLN1B_MYCTU'
# a68<-desc.reordered=='sp_P99063_ODPB_STAAN'
# a69<-desc.reordered=='tr_Q9CHQ7_Q9CHQ7_LACLA'
# a70<-desc.reordered=='sp_P0ABK2_CYDB_ECOLI'
# a71<-desc.reordered=='sp_P32166_MENA_ECOLI'
# a72<-desc.reordered=='DEG10180287'
# a73<-desc.reordered=='sp_P0C1R8_MRAY_STAAU'
# a74<-desc.reordered=='sp_Q56584_NQRF_VIBAL'
# a75<-desc.reordered=='DEG10180599'
# a76<-desc.reordered=='sp_P00804_LSPA_ECOLI'
# a77<-desc.reordered=='sp_P60089_ODPA_STAAM'
# a78<-desc.reordered=='sp_Q9RDT5_WALR_STAAU'
# a79<-desc.reordered=='DEG10180117'
# a80<-desc.reordered=='sp_P0A6Z3_HTPG_ECOLI'
# a81<-desc.reordered=='DEG10180529'
# a82<-desc.reordered=='DEG10180054'
# a83<-desc.reordered=='sp_P0AG30_RHO_ECOLI'
# a84<-desc.reordered=='DEG10180586'
# a85<-desc.reordered=='sp_P0AEK2_FABG_ECOLI'
# a86<-desc.reordered=='tr_A0A3G5FN01_A0A3G5FN01_ECOLX'
# a87<-desc.reordered=='sp_P0ABU0_MENB_ECOLI'
# a88<-desc.reordered=='sp_P0A9M0_LON_ECOLI'
# a89<-desc.reordered=='sp_P0AGA2_SECY_ECOLI'
# a90<-desc.reordered=='DEG10180342'
# a91<-desc.reordered=='sp_P11880_MURF_ECOLI'
# a92<-desc.reordered=='DEG10180349'
# a93<-desc.reordered=='DEG10180510'
# a94<-desc.reordered=='sp_P0ADG7_IMDH_ECOLI'
# myCols[a1]='#000000'
# myCols[a2]='#012C58'
# myCols[a3]='#1CE6FF'
# myCols[a4]='#FF34FF'
# myCols[a5]='#FF4A46'
# myCols[a6]='#008941'
# myCols[a7]='#006FA6'
# myCols[a8]='#A30059'
# myCols[a9]='#FFDBE5'
# myCols[a10]='#7A4900'
# myCols[a11]='#0000A6'
# myCols[a12]='#63FFAC'
# myCols[a13]='#B79762'
# myCols[a14]='#004D43'
# myCols[a15]='#8FB0FF'
# myCols[a16]='#997D87'
# myCols[a17]='#5A0007'
# myCols[a18]='#809693'
# myCols[a19]='#7ED379'
# myCols[a20]='#1B4400'
# myCols[a21]='#4FC601'
# myCols[a22]='#372101'
# myCols[a23]='#4A3B53'
# myCols[a24]='#FF2F80'
# myCols[a25]='#11915A'
# myCols[a26]='#BA0900'
# myCols[a27]='#6B7900'
# myCols[a28]='#00C2A0'
# myCols[a29]='#FFAA92'
# myCols[a30]='#FF90C9'
# myCols[a31]='#B903AA'
# myCols[a32]='#D16100'
# myCols[a33]='#DDEFFF'
# myCols[a34]='#000035'
# myCols[a35]='#7B4F4B'
# myCols[a36]='#A1C299'
# myCols[a37]='#300018'
# myCols[a38]='#0AA6D8'
# myCols[a39]='#013349'
# myCols[a40]='#00846F'
# myCols[a41]='#3B5DFF'
# myCols[a42]='#FFB500'
# myCols[a43]='#C2FFED'
# myCols[a44]='#A079BF'
# myCols[a45]='#CC0744'
# myCols[a46]='#C0B9B2'
# myCols[a47]='#C2FF99'
# myCols[a48]='#001E09'
# myCols[a49]='#00489C'
# myCols[a50]='#6F0062'
# myCols[a51]='#0CBD66'
# myCols[a52]='#EEC3FF'
# myCols[a53]='#456D75'
# myCols[a54]='#B77B68'
# myCols[a55]='#7A87A1'
# myCols[a56]='#788D66'
# myCols[a57]='#885578'
# myCols[a58]='#FAD09F'
# myCols[a59]='#FF8A9A'
# myCols[a60]='#D157A0'
# myCols[a61]='#BEC459'
# myCols[a62]='#456648'
# myCols[a63]='#0086ED'
# myCols[a64]='#886F4C'
# myCols[a65]='#34362D'
# myCols[a66]='#B4A8BD'
# myCols[a67]='#00A6AA'
# myCols[a68]='#452C2C'
# myCols[a69]='#636375'
# myCols[a70]='#A3C8C9'
# myCols[a71]='#FF913F'
# myCols[a72]='#938A81'
# myCols[a73]='#575329'
# myCols[a74]='#00FECF'
# myCols[a75]='#B05B6F'
# myCols[a76]='#8CD0FF'
# myCols[a77]='#3B9700'
# myCols[a78]='#04F757'
# myCols[a79]='#C8A1A1'
# myCols[a80]='#1E6E00'
# myCols[a81]='#7900D7'
# myCols[a82]='#A77500'
# myCols[a83]='#6367A9'
# myCols[a84]='#A05837'
# myCols[a85]='#6B002C'
# myCols[a86]='#772600'
# myCols[a87]='#D790FF'
# myCols[a88]='#9B9700'
# myCols[a89]='#549E79'
# myCols[a90]='#FFF69F'
# myCols[a91]='#201625'
# myCols[a92]='#72418F'
# myCols[a93]='#BC23FF'
# myCols[a94]='#99ADC0'


# for 616 targets
a1<-desc.reordered=='DEG10180276'
myCols[a1]='#012C58'
a2<-desc.reordered=='DEG10180273'
myCols[a2]='#1CE6FF'
a3<-desc.reordered=='DEG10180377'
myCols[a3]='#FF34FF'
a4<-desc.reordered=='DEG10180080'
myCols[a4]='#FF4A46'
a5<-desc.reordered=='DEG10180475'
myCols[a5]='#008941'
a6<-desc.reordered=='DEG10180498'
myCols[a6]='#006FA6'
a7<-desc.reordered=='DEG10180490'
myCols[a7]='#A30059'
a8<-desc.reordered=='DEG10180016'
myCols[a8]='#FFDBE5'
a9<-desc.reordered=='DEG10180400'
myCols[a9]='#7A4900'
a10<-desc.reordered=='DEG10180012'
myCols[a10]='#0000A6'
a11<-desc.reordered=='DEG10180013'
myCols[a11]='#63FFAC'
a12<-desc.reordered=='DEG10180155'
myCols[a12]='#B79762'
a13<-desc.reordered=='DEG10180156'
myCols[a13]='#004D43'
a14<-desc.reordered=='DEG10180150'
myCols[a14]='#8FB0FF'
a15<-desc.reordered=='DEG10180335'
myCols[a15]='#997D87'
a16<-desc.reordered=='DEG10180334'
myCols[a16]='#5A0007'
a17<-desc.reordered=='DEG10180204'
myCols[a17]='#809693'
a18<-desc.reordered=='DEG10180200'
myCols[a18]='#7ED379'
a19<-desc.reordered=='DEG10180208'
myCols[a19]='#1B4400'
a20<-desc.reordered=='DEG10180595'
myCols[a20]='#4FC601'
a21<-desc.reordered=='DEG10180597'
myCols[a21]='#3B5DFF'
a22<-desc.reordered=='DEG10180596'
myCols[a22]='#4A3B53'
a23<-desc.reordered=='DEG10180590'
myCols[a23]='#FF2F80'
a24<-desc.reordered=='DEG10180110'
myCols[a24]='#11915A'
a25<-desc.reordered=='DEG10180113'
myCols[a25]='#BA0900'
a26<-desc.reordered=='DEG10180443'
myCols[a26]='#6B7900'
a27<-desc.reordered=='DEG10180522'
myCols[a27]='#00C2A0'
a28<-desc.reordered=='DEG10180296'
myCols[a28]='#FFAA92'
a29<-desc.reordered=='DEG10180529'
myCols[a29]='#FF90C9'
a30<-desc.reordered=='DEG10180069'
myCols[a30]='#B903AA'
a31<-desc.reordered=='DEG10180066'
myCols[a31]='#D16100'
a32<-desc.reordered=='DEG10180062'
myCols[a32]='#DDEFFF'
a33<-desc.reordered=='DEG10180115'
myCols[a33]='#000035'
a34<-desc.reordered=='DEG10180114'
myCols[a34]='#7B4F4B'
a35<-desc.reordered=='DEG10180524'
myCols[a35]='#A1C299'
a36<-desc.reordered=='DEG10180241'
myCols[a36]='#300018'
a37<-desc.reordered=='DEG10180240'
myCols[a37]='#0AA6D8'
a38<-desc.reordered=='DEG10180242'
myCols[a38]='#013349'
a39<-desc.reordered=='DEG10180244'
myCols[a39]='#00846F'
a40<-desc.reordered=='DEG10180246'
myCols[a40]='#372101'
a41<-desc.reordered=='DEG10180360'
myCols[a41]='#FFB500'
a42<-desc.reordered=='DEG10180361'
myCols[a42]='#C2FFED'
a43<-desc.reordered=='DEG10180366'
myCols[a43]='#A079BF'
a44<-desc.reordered=='DEG10180367'
myCols[a44]='#CC0744'
a45<-desc.reordered=='DEG10180098'
myCols[a45]='#C0B9B2'
a46<-desc.reordered=='DEG10180099'
myCols[a46]='#C2FF99'
a47<-desc.reordered=='DEG10180568'
myCols[a47]='#001E09'
a48<-desc.reordered=='DEG10180564'
myCols[a48]='#00489C'
a49<-desc.reordered=='DEG10180397'
myCols[a49]='#6F0062'
a50<-desc.reordered=='DEG10180396'
myCols[a50]='#0CBD66'
a51<-desc.reordered=='DEG10180393'
myCols[a51]='#EEC3FF'
a52<-desc.reordered=='DEG10180392'
myCols[a52]='#456D75'
a53<-desc.reordered=='DEG10180391'
myCols[a53]='#B77B68'
a54<-desc.reordered=='DEG10180029'
myCols[a54]='#7A87A1'
a55<-desc.reordered=='DEG10180021'
myCols[a55]='#788D66'
a56<-desc.reordered=='DEG10180141'
myCols[a56]='#885578'
a57<-desc.reordered=='DEG10180145'
myCols[a57]='#FAD09F'
a58<-desc.reordered=='DEG10180418'
myCols[a58]='#FF8A9A'
a59<-desc.reordered=='DEG10180148'
myCols[a59]='#D157A0'
a60<-desc.reordered=='DEG10180325'
myCols[a60]='#BEC459'
a61<-desc.reordered=='DEG10180323'
myCols[a61]='#456648'
a62<-desc.reordered=='DEG10180321'
myCols[a62]='#0086ED'
a63<-desc.reordered=='DEG10180318'
myCols[a63]='#886F4C'
a64<-desc.reordered=='DEG10180190'
myCols[a64]='#34362D'
a65<-desc.reordered=='DEG10180218'
myCols[a65]='#B4A8BD'
a66<-desc.reordered=='DEG10180219'
myCols[a66]='#00A6AA'
a67<-desc.reordered=='DEG10180197'
myCols[a67]='#452C2C'
a68<-desc.reordered=='DEG10180198'
myCols[a68]='#636375'
a69<-desc.reordered=='DEG10180210'
myCols[a69]='#A3C8C9'
a70<-desc.reordered=='DEG10180214'
myCols[a70]='#FF913F'
a71<-desc.reordered=='DEG10180352'
myCols[a71]='#938A81'
a72<-desc.reordered=='DEG10180357'
myCols[a72]='#575329'
a73<-desc.reordered=='DEG10180356'
myCols[a73]='#00FECF'
a74<-desc.reordered=='DEG10180608'
myCols[a74]='#B05B6F'
a75<-desc.reordered=='DEG10180359'
myCols[a75]='#8CD0FF'
a76<-desc.reordered=='DEG10180580'
myCols[a76]='#3B9700'
a77<-desc.reordered=='DEG10180581'
myCols[a77]='#04F757'
a78<-desc.reordered=='DEG10180584'
myCols[a78]='#C8A1A1'
a79<-desc.reordered=='DEG10180585'
myCols[a79]='#1E6E00'
a80<-desc.reordered=='DEG10180589'
myCols[a80]='#7900D7'
a81<-desc.reordered=='DEG10180107'
myCols[a81]='#A77500'
a82<-desc.reordered=='DEG10180452'
myCols[a82]='#6367A9'
a83<-desc.reordered=='DEG10180513'
myCols[a83]='#A05837'
a84<-desc.reordered=='P0AAI3'
myCols[a84]='#6B002C'
a85<-desc.reordered=='DEG10180072'
myCols[a85]='#772600'
a86<-desc.reordered=='DEG10180019'
myCols[a86]='#D790FF'
a87<-desc.reordered=='DEG10180255'
myCols[a87]='#9B9700'
a88<-desc.reordered=='DEG10180253'
myCols[a88]='#549E79'
a89<-desc.reordered=='DEG10180315'
myCols[a89]='#FFF69F'
a90<-desc.reordered=='DEG10180313'
myCols[a90]='#201625'
a91<-desc.reordered=='DEG10180312'
myCols[a91]='#72418F'
a92<-desc.reordered=='DEG10180311'
myCols[a92]='#BC23FF'
a93<-desc.reordered=='DEG10180533'
myCols[a93]='#99ADC0'
a94<-desc.reordered=='DEG10180537'
myCols[a94]='#3A2465'
a95<-desc.reordered=='DEG10180122'
myCols[a95]='#922329'
a96<-desc.reordered=='DEG10180509'
myCols[a96]='#5B4534'
a97<-desc.reordered=='DEG10180539'
myCols[a97]='#FDE8DC'
a98<-desc.reordered=='DEG10180222'
myCols[a98]='#404E55'
a99<-desc.reordered=='DEG10180552'
myCols[a99]='#0089A3'
a100<-desc.reordered=='DEG10180554'
myCols[a100]='#CB7E98'
a101<-desc.reordered=='DEG10180281'
myCols[a101]='#A4E804'
a102<-desc.reordered=='DEG10180382'
myCols[a102]='#324E72'
a103<-desc.reordered=='DEG10180388'
myCols[a103]='#6A3A4C'
a104<-desc.reordered=='DEG10180389'
myCols[a104]='#83AB58'
a105<-desc.reordered=='DEG10180038'
myCols[a105]='#001C1E'
a106<-desc.reordered=='DEG10180034'
myCols[a106]='#D1F7CE'
a107<-desc.reordered=='DEG10180133'
myCols[a107]='#004B28'
a108<-desc.reordered=='DEG10180132'
myCols[a108]='#C8D0F6'
a109<-desc.reordered=='DEG10180220'
myCols[a109]='#A3A489'
a110<-desc.reordered=='DEG10180042'
myCols[a110]='#806C66'
a111<-desc.reordered=='DEG10180047'
myCols[a111]='#222800'
a112<-desc.reordered=='DEG10180044'
myCols[a112]='#BF5650'
a113<-desc.reordered=='DEG10180185'
myCols[a113]='#E83000'
a114<-desc.reordered=='DEG10180180'
myCols[a114]='#66796D'
a115<-desc.reordered=='DEG10180267'
myCols[a115]='#DA007C'
a116<-desc.reordered=='DEG10180264'
myCols[a116]='#FF1A59'
a117<-desc.reordered=='DEG10180342'
myCols[a117]='#8ADBB4'
a118<-desc.reordered=='DEG10180344'
myCols[a118]='#1E0200'
a119<-desc.reordered=='DEG10180363'
myCols[a119]='#5B4E51'
a120<-desc.reordered=='DEG10180500'
myCols[a120]='#C895C5'
a121<-desc.reordered=='DEG10180488'
myCols[a121]='#320033'
a122<-desc.reordered=='DEG10180483'
myCols[a122]='#FF6832'
a123<-desc.reordered=='DEG10180482'
myCols[a123]='#66E1D3'
a124<-desc.reordered=='DEG10180486'
myCols[a124]='#CFCDAC'
a125<-desc.reordered=='DEG10180004'
myCols[a125]='#D0AC94'
b0<-desc.reordered=='DEG10180002'
myCols[b0]='#000000'
b1<-desc.reordered=='DEG10180164'
myCols[b1]='#012C58'
b2<-desc.reordered=='O86737'
myCols[b2]='#1CE6FF'
b3<-desc.reordered=='DEG10180304'
myCols[b3]='#FF34FF'
b4<-desc.reordered=='DEG10180300'
myCols[b4]='#FF4A46'
b5<-desc.reordered=='DEG10180301'
myCols[b5]='#008941'
b6<-desc.reordered=='DEG10180302'
myCols[b6]='#006FA6'
b7<-desc.reordered=='DEG10180303'
myCols[b7]='#A30059'
b8<-desc.reordered=='DEG10180008'
myCols[b8]='#FFDBE5'
b9<-desc.reordered=='DEG10180545'
myCols[b9]='#7A4900'
b10<-desc.reordered=='DEG10180236'
myCols[b10]='#0000A6'
b11<-desc.reordered=='DEG10180238'
myCols[b11]='#63FFAC'
b12<-desc.reordered=='DEG10180135'
myCols[b12]='#B79762'
b13<-desc.reordered=='DEG10180549'
myCols[b13]='#004D43'
b14<-desc.reordered=='DEG10180438'
myCols[b14]='#8FB0FF'
b15<-desc.reordered=='DEG10180439'
myCols[b15]='#997D87'
b16<-desc.reordered=='DEG10180435'
myCols[b16]='#5A0007'
b17<-desc.reordered=='DEG10180125'
myCols[b17]='#809693'
b18<-desc.reordered=='DEG10180289'
myCols[b18]='#7ED379'
b19<-desc.reordered=='DEG10180285'
myCols[b19]='#1B4400'
b20<-desc.reordered=='DEG10180286'
myCols[b20]='#4FC601'
b21<-desc.reordered=='DEG10180283'
myCols[b21]='#3B5DFF'
b22<-desc.reordered=='DEG10180282'
myCols[b22]='#4A3B53'
b23<-desc.reordered=='DEG10180050'
myCols[b23]='#FF2F80'
b24<-desc.reordered=='DEG10180051'
myCols[b24]='#11915A'
b25<-desc.reordered=='DEG10180052'
myCols[b25]='#BA0900'
b26<-desc.reordered=='DEG10180054'
myCols[b26]='#6B7900'

myBG <- myCols

##############################################################
# Read the 2f file of fastaID and description
dir1 <- "/Users/gvandova/Dropbox/Computational_projects/TargetMiningGenomes/Phylogeny/"
filename1 <- "KS.14.10kb.fasta.filtered.phyla" # to get colors by phyla for KSs
# filename1 <- "KS.616.10kb.fasta.filtered.phyla" # to get colors by phyla for KSs
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

mywidth=4; myheight=6 #for 14 targets, pos set
# mywidth=6; myheight=8 #for 14 targets, show labels
# mywidth=6; myheight=30 #for 119 targets rooted tree
# mywidth=20; myheight=86 #for 609 targets rooted tree

outfile <- paste(dir, filename, ".refset.png", sep="")
title <- "KS.14.10kb"
pdf(file=outfile, width=mywidth, height=myheight)
plot(MyTree.ladderized, main=title, font=.1, type="phylogram", edge.color="gray",
     edge.width=.5, show.tip.label=F, open.angle=5) # for rooted tree
# plot(MyTree, font=1, type="unrooted", edge.color="gray", edge.width=.5, show.tip.label=F, open.angle=5) # for unrooted tree
tiplabels(pch=21, cex=.5, col=myCols, bg=myBG)# pch=21 circles, colored by target, refset
# tiplabels(pch=21, cex=.3, col=myCols, bg=myBG)# pch=21 circles, colored by target, all other
# tiplabels(pch=21, cex=.5, col=myCols1, bg=myBG1)# colored by phyla

# Print labels
# tiplabels(MyTree$tip.label, cex=.2, frame="none", adj=0) # full label
# tiplabels(phyla.reordered, cex=.1, frame="none", adj=0) # phyla label
# tiplabels(phyla.reordered, cex=.5, frame="none", adj=0) # phyla  label 609

# tiplabels(desc2.reordered, cex=0.2, frame="none", adj=0) # short descr label
tiplabels(desc3.reordered, cex=0.3, frame="none", adj=0) # to highlight mibig ref sequences

# nodelabels(MyTree$node.label, frame="none", cex=.1)
# edgelabels(MyTree$edge.label, frame="none", cex=.2)

add.scale.bar(cex=.5, lwd=.5)
dev.off()
