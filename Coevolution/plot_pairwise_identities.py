#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pyplot
import sys
import pylab
from matplotlib.pyplot import *  # This is for the legend to work
from collections import defaultdict

ks_pair_identities = defaultdict(list)
target_pair_identities = defaultdict(list)
pairs = []

# pairwis_filename = "pairwise_identities.mibig.out"
pairwise_filename = "pairwise_identities.92.5kb.out"
with open(pairwise_filename, 'r') as f:
    data = f.readlines()
    data = map(lambda x: x.strip(), data)
    for line in data:
        print line
        gene1, gene2, len1, len2, ks_pair_identity, target_pair_identity, d = line.split("\t")
        target = gene1.split("|")[1]
        print gene1, gene2, target
        ks_pair_identities[target].append(float(ks_pair_identity))
        target_pair_identities[target].append(float(target_pair_identity))
        pairs.append((gene1,gene2))

# ks_pair_identities = map(float, ks_pair_identities)
# target_pair_identities = map(float, target_pair_identities)

print len(ks_pair_identities), len(target_pair_identities)

colors = [
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF",
    "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC",
    "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#FEFFE6", "#1B4400",
    "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0",
    "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
    "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF",
    "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF",
    "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
    "#BEC459", "#456648", "#0086ED", "#886F4C",

    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C",
    "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF",
    "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837",
    "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F",
    "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
    "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
    "#83AB58", "#001C1E", "#D1F7CE", "#004B28",
    "#C8D0F6", "#A3A489", "#806C66", "#222800",
    "#BF5650", "#E83000", "#66796D", "#DA007C",
    "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
    "#C895C5", "#320033", "#FF6832", "#66E1D3",
    "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
]

targets = [
    'sp_P17109_MEND_ECOLI', 'sp_P0A6K3_DEF_ECOLI', 'mfR6_squalene_synthase',
    'sp_P16659_SYP_ECOLI', 'sp_P05041_PABB_ECOLI', 'sp_P0A8M3_SYT_ECOLI',
    'sp_P0A6W3_MRAY_ECOLI', 'sp_P0A7Z4_RPOA_ECOLI', 'sp_P60785_LEPA_ECOLI',
    'sp_P10443_DPO3A_ECOLI', 'sp_Q9RDT5_WALR_STAAU', 'sp_P21889_SYD_ECOLI',
    'sp_P0A6G7_CLPP_ECOLI', 'tr_E2QIY7_E2QIY7_ECOLX', 'sp_P0A6M8_EFG_ECOLI',
    'sp_P77781_MENI_ECOLI', 'sp_P0A8L1_SYS_ECOLI', 'sp_P04079_GUAA_ECOLI',
    'sp_P07862_DDLB_ECOLI', 'mpaF_IMDH', 'sp_P0CE47_EFTU1_ECOLI',
    'sp_P0AGA2_SECY_ECOLI', 'sp_P28305_PABC_ECOLI', 'sp_P0A725_LPXC_ECOLI',
    'sp_P10408_SECA_ECOLI', 'sp_P0AG30_RHO_ECOLI', 'sp_P0AGB6_RPOE_ECOLI',
    'sp_P08312_SYFA_ECOLI', 'sp_P0AEK2_FABG_ECOLI', 'sp_Q9RDT3_WALK_STAAU',
    'sp_P0A749_MURA_ECOLI', 'sp_P45568_DXR_ECOLI', 'sp_P0A9M0_LON_ECOLI',
    'sp_P08373_MURB_ECOLI', 'sp_P0A6B4_ALR1_ECOLI', 'beta_lactamase',
    'sp_P17169_GLMS_ECOLI', 'sp_P06986_HIS8_ECOLI', 'sp_P0AGJ9_SYY_ECOLI',
    'sp_P0AAI3_FTSH_ECOLI', 'sp_Q75R59_NQRF_VIBAN', 'sp_P0A6N4_EFP_ECOLI',
    'sp_P07003_POXB_ECOLI', 'sp_P00903_PABA_ECOLI', 'sp_P32166_MENA_ECOLI',
    'sp_P0ABU0_MENB_ECOLI', 'sp_P0C0V0_DEGP_ECOLI', 'sp_P37353_MENE_ECOLI',
    'sp_P03007_DPO3E_ECOLI', 'sp_P18335_ARGD_ECOLI'
]

target_to_color = dict([(k, v) for k, v in zip(targets, colors)])
# target_to_color = {
#         "sp_P0A6K3_DEF_ECOLI": '#CFCDAC',
#         "SalI_beta_proteasome": 'lightblue',
#         "GriR_DnaN": "cyan",
#         "EF-Tu": "midnightblue",
#         "PtmP3_FabB-F": "r",
#         "BatG_FabI": "darkmagenta",
#         "GyrB-R": "magenta",
#         "mupM_Ile-tRNA-syn": "brown",
#         "borI_Thr-tRNA-syn": "lightpink",
#         "agnB2_Leu-tRNA-syn": "orange",
#         "rubR1_TIF": "green",
#         "Ind0_Trp-tRNA-syn": "lightgreen"}

target_to_name = dict([(k, k) for k in targets])
for target in ks_pair_identities:
    # plt.scatter(ks_pair_identities, target_pair_identities, color='r', s=1)
    if target not in target_to_color:
        color = "black"
    else:
        color = target_to_color[target]
    if target not in target_to_name:
        label=None
    else:
        label = target_to_name[target]
    plt.scatter(ks_pair_identities[target], target_pair_identities[target],
                color=color, s=10, label=label)
    plt.xlabel('KS1-KS2 identity', size=10)
    plt.ylabel('Target1-Target2 identity', size=10)
    plt.legend(loc=4, fontsize="xx-small")


    # y = [2.56422, 3.77284, 3.52623, 3.51468, 3.02199]
    # z = [0.15, 0.3, 0.45, 0.6, 0.75]
    # n = [58, 651, 393, 203, 123]
    #
    # fig, ax = plt.subplots()

# #Doesn't work
# for i, txt in enumerate(pairs):
#     plt.annotate(txt, (ks_pair_identities[i], target_pair_identities[i]))


plt.xlim([0, 100])
plt.ylim([0, 100])
# plt.savefig('mibig.png')
plt.savefig('92.5kb.png', dpi=400)
