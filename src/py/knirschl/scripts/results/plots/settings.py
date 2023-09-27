from constants import *
import matplotlib.pyplot as plt

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

COLORS = {
    GENERAX: "#e6194B",
    RAXMLNG: "#800000",
    FASTME: "#f58231",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, APRO)]: "#f58231", #"#000075",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, MAD)]: "#469990", #"#469990",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, ALL)]: "#000075",
    SPEARFISH_MAP[BA_FASTME]: "#4363d8",
    APRO: "#000075",
    MAD: "#469990",
    ALL: "#4363d8",
    OTHER: "#542788",
    "Spearfish (all": "#4363d8",
    None: None
}

MARKERS = {
    GENERAX: "o",
    RAXMLNG: "x",
    FASTME: "d",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, APRO)]: "o",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, MAD)]: "x",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, ALL)]: "d",
    SPEARFISH_MAP[BA_FASTME]: "d",
    APRO: "o",
    MAD: "x",
    ALL: "d",
    OTHER: "o",
    "Spearfish (all": "d",
    None: None
}

MARKERSIZE = 2.5