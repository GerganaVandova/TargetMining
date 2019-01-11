MIBiG Target-In-Cluster Analysis Results
========================================
The experiment went as follows:

1. Download json and fasta files from MIBiG http://mibig.secondarymetabolites.org/download.html
2. Parse the MIBiG json files to find putative targets for each MIBiG entry ("Molecular target" on the website).
   These are free text and have typos etc.
3. Search uniprot for a set of proteins matching this description. Sometimes it's a single gene, sometimes it's broad, like "cell wall".
   My general threshold was to return no more than 500 proteins.
   Sometimes I use uniref to achieve this (it clusters by identity and returns only proteins with < 50% identity).
4. Make blast dbs from each of the proteins in the clusters in the MIBiG fasta file (`MIBiG_prot_seqs_1.3.fasta`)
   `for f in BGC*fasta; do makeblastdb -in $f -dbtype prot -parse_seqids -out $(basename $f .fasta); done`
5. For each protein in my set of uniprot fastas, check if there is a homolog in the cluster with blastp.
   `python run_blast.py` (which uses `blastp -query uniprot_fasta/{fasta} -db blastdb/{cluster} -out out/{cluster}.out -evalue 1e-50 -outfmt 6`)

NOT FOUND
BGC0000023
----------
Target: "cholesterol biosynthesis"
Query: sp|Q9Y7D5|LOVF_ASPTE
Hit: BGC0000023|c1|8559-14162|-|no_locus_tag|polyketide_synthase_AufC|CAO98847.1
evalue: 1e-140
False positive. LOVF is "Lovastatin diketide synthase" (http://www.uniprot.org/uniprot/Q9Y7D5)
HMGR was in the set of uniprot genes but no homolog was found.


FOUND
BGC0000031
----------
Target: "Threonyl-tRNA synthetase"
Query: SYTC_HUMAN
Hit: BGC0000031|c1|57939-59966|+|no_locus_tag|Threonyl-tRNA_synthetase|CAE45679.1
evalue: 0
Pretty obvious since it's annotated correctly in MIBiG.


FOUND
BGC0000104
----------
Target: "inosine monophosphate dehydrogenase"
Query: IMDH1_HUMAN
Hit: BGC0000104|c1|18958-20657|-|no_locus_tag|putative_inosine_monophosphate_dehydrogenase|ADY00133.1
evalue: 0
One of our true positives.


BGC0000168
----------
Target: "undecaprenyl diphosphate synthase"
Query: sp|D7PHZ2|VRTA_PENAE
Hit: BGC0000168|c1|20532-26109|+|no_locus_tag|VrtA|ADI24926.1
evalue: 0
False positive. The query genes are just the genes from the biosynthetic cluster.


Found
BGC0000182
----------
Target: "Isoleucyl tRNA synthetase"
Query: SYIC_HUMAN
Hit: BGC0000182|c1|54444-57536|+|no_locus_tag|MupM|AAM12927.2
evalue: 0
MupM is isoleucyl tRNA synthase (https://en.wikipedia.org/wiki/Mupirocin)

FOUND
BGC0000345
----------
Target: "proeasome" [proteasome]
Query: UniRef50_Q53079
Hit: BGC0000345|c1|2040-2897|+|no_locus_tag|EpnC|AHB38505.1
evalue: 1e-100
This is plausible.
The hit has proteasome domains (https://www.ncbi.nlm.nih.gov/protein/563322452/)


BGC0000447
----------
Target: "cell membrane"
Query: UniRef50_O31101	UniRef50_P58411
Hit: BGC0000447|c1|65838-66989|+|no_locus_tag|macrolide-specific_efflux_protein_macA|CCJ67641.1	BGC0000447|c1|1-1440|-|no_locus_tag|multidrug/solvent_efflux_pump_outer_membrane_protein_mepC|CCJ67633.1
evalue: 1e-90	0
This was one of the extremely broad targets groups.


Found
BGC0000615
----------
Target: "ribosome"
Query: UniRef50_A2R994
    Hit: BGC0000615|c1|31692-33803|+|no_locus_tag|Tpa10|ADO67785.1
evalue: 1e-100
This one seems plausible.
The NCBI entry (https://www.ncbi.nlm.nih.gov/protein/ADO67785.1) is a ribosomal elongation factor (bacterial).
The entry was from a paper: "Isolation and Characterization of the Gene Cluster for Biosynthesis of the Thiopeptide Antibiotic TP-1161"


FOUND
BGC0000832
----------
Target: "GyrB, subunit B protein of DNA gyrase"
Query: GYRB_ECOLI
Hit: BGC0000832|c1|31227-32944|+|no_locus_tag|clorobiocin-resistant_gyrase_B|AAN65247.1
evalue: 0
This one is obvious.

FOUND
BGC0000833
----------
Target: "GyrB, subunit B protein of DNA gyrase"
Query: GYRB_ECOLI
Hit: BGC0000833|c1|33611-35644|+|no_locus_tag|DNA_gyrase_subunit_B|AAO47225.2
evalue: 0
This one is obvious.

FOUND
BGC0000834
----------
Target: "GyrB, subunit B protein of DNA gyrase"
Query: GYRB_ECOLI
Hit: BGC0000834|c1|21185-23218|+|no_locus_tag|GyrB_R|AFI47646.1
evalue: 0
This one is obvious


BGC0000849
----------
Target: "SCO6265"
Query: tr|Q7AKF1|Q7AKF1_STRCO
Hit: BGC0000849|c1|1-648|-|no_locus_tag|gamma-butyrolactone_binding_protein|CAB60184.1
evalue: 1e-160
SCO6265 is Gamma-butyrolactone binding protein (http://www.uniprot.org/uniprot/Q7AKF1)


BGC0000956
----------
Target: "beta-subunit of acetyl-CoA carboxylase"
Query: UniRef50_A2CAG3
Hit: BGC0000956|c1|24422-25336|+|no_locus_tag|AdmT|AAO39114.1
evalue: 0
The protein is "acetyl-coenzyme A carboxylase carboxyl transferase subunit beta" (https://www.ncbi.nlm.nih.gov/protein/AAO39114.1)


BGC0000971
----------
Target: "proteasome inhibitor"
Query: UniRef50_Q53079
Hit: BGC0000971|c1|14472-15317|-|no_locus_tag|20S_proteasome_beta-subunit|CBW54670.1
evalue: 1e-100
Another very broad category, but possible.


BGC0001041
----------
Target: "20S proteasome"
Query: sp|P9WHT9|PSB_MYCTU
Hit: BGC0001041|c1|3959-4807|+|Strop_1015|20S_proteasome,_A_and_B_subunits|ABP53490.1
evalue: 1e-100
This seems plausible.


BGC0001067
----------
Target: "methionine aminopeptidase 2"
Query: sp|P50579|MAP2_HUMAN
Hit: BGC0001067|c1|11015-12951|-|AFUA_8G00410|methionine_aminopeptidase,_type_II,_putative|EAL85125.1
evalue: 1e-100
This is one of our true positives.


BGC0001082
----------
Target: "BK channels"
Query: sp|Q9C446|PAXG_PENPX
Hit: BGC0001082|c1|1-1294|-|no_locus_tag|PaxG|AAK11531.1
evalue: 0
This is a false positive due to pulling in genes labeled "potassium channel blocker"
(they biosynthesize paxilline, this potassium channel blocker). (http://www.uniprot.org/uniprot/Q9C446)


BGC0001099
----------
Target: "FabI"
Query: sp|Q2P9J6|FABV_XANOM
Hit: BGC0001099|c1|52666-53856|+|no_locus_tag|BatG|ADD82948.1
evalue: 1e-150
The NCBI entry for BatG is from the paper "Isolation and purification of a new kalimantacin/batumin-related
polyketide antibiotic and elucidation of its biosynthesis gene cluster"
(https://www.ncbi.nlm.nih.gov/protein/ADD82948.1)
FabI is Enoyl-[acyl-carrier-protein] reductase (NADH) (https://www.ncbi.nlm.nih.gov/protein/ABS73537.1)
BatG is trans-2-enoyl-CoA reductase.
BatG is about 22% identical to FabI, so not that close (58/258).
Still, this is very plausible.

```
CLUSTAL multiple sequence alignment by MUSCLE (3.8)

FabI            --MNFSLEG----------------------RNIVVMGVANKRSIAWGIARSLHEAGARL
BatG            MIVNPRVKGFICTTAHPAGCRANVDEQIRFIREQAPIANAPKRVLVIGASTG-YGLASRI
                  :*  ::*                      *: . :. * ** :. * : . :  .:*:

FabI            IFTYA------GERLEKSV-HDLAATLERNDSIILPC--------------DVTNDAEIE
BatG            TAAFGCNARTIGVFFEKPASRNRTASAGWYNSAAFQCAADNAGLYAKSINGDAFSDAVKE
                  ::.      *  :**.. .: :*:  . :*  : *              *. .**  *

FabI            ACFASIKEQVGVIHGIAHCIAFANKEELVGEYLNTN----------REGFLLAHN-ISSY
BatG            KTIEMIKADLGQIDLLVYSLASPRRLHPVTGTLHSSVLKPIGKTVTQIGLDTDRELIKSF
                  :  ** ::* *  :.:.:* ...   *   *::.          . *:   .: *.*:

FabI            SLTAVAK----------------------AARPIMTEGGSIVTLTYLGGERVVSNY--NV
BatG            TLQPAVQQEIDDTVVVMGADDWERWVHQLSEAGVLAPGCKTTVYTYIGEKVTRDIYWDGT
                :* ...:                      :   ::: * . .. **:* : . . *  ..

FabI            MGVAKASLDASVRYLAADLGKENIRV------NSISAGPIRTL-------------SAKG
BatG            IGAAKKDLDRAAAMLSSNGVEASVSVLKAVVTQSSAAIPVMPLYLALLFKLMKQDGSHEG
                :*.** .** :.  *:::  : .: *      :* :* *: .*             * :*

FabI            ISD----------FNS--------------------ILKEIEERAPLRRTTTPEEVGDTA
BatG            CIEQIYRLFSECLYNTDPRLDEGGRNRVDDRETRHEIQAEVEKLWPQVTTENLNSISDFQ
                  :          :*:                    *  *:*:  *   * . :.:.*  

FabI            AF---LFSDLSRGITGENLHVDSGFHITAR---------
BatG            GFRAEFLKLFGFGLT--EIDYETDFDVDVKINGLLDLTQ
                .*   ::. :. *:*  ::  ::.* : ..       
```

BGC0001140
----------
Target: "Fatty acid synthases FabF & FabH"
Query: UniRef50_O34340
Hit: BGC0001140|c1|37763-38977|+|no_locus_tag|PtmP3|ACS13710.1
evalue: 1e-80
The hit is has a FabF region (https://www.ncbi.nlm.nih.gov/protein/ACS13710.1)
"PtmP3/PtnP3 and FabF confer PTM and PTN resistance by target replacement and target modification, respectively."
(http://www.genengnews.com/gen-news-highlights/bacterial-resistance-mechanism-exposed-new-antibiotics-may-follow/81249538)


BGC0001155
----------
Target: "EF-Tu"
Query: sp|P0CE47|EFTU1_ECOLI
Hit: BGC0001155|c1|25094-26287|+|no_locus_tag|EF_Tu|AGY49600.1
evalue: 0
This is obvious.


BGC0001156
----------
Target: "FabF, FabH"
Query: UniRef50_O34340
Hit: BGC0001156|c1|32358-33572|+|no_locus_tag|PtnP3|ADD83010.1
evalue: 1e-80
See BGC0001140.


BGC0001339
----------
Target: "squalene synthase"
Query: sp|P37268|FDFT_HUMAN
Hit: BGC0001339|c1|66436-67748|+|no_locus_tag|squalene_synthase|AMY15074.1
FDFT1 is squalene synthase (http://www.uniprot.org/uniprot/P37268).
