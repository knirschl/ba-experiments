from constants import *

COLORS = {
    GENERAX: "#e6194B",
    RAXML: "#800000",
    FASTME: "#f58231",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, APRO)]: "#000075",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, MAD)]: "#469990",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, ALL)]: "#4363d8",
    SPEARFISH_MAP[BA_FASTME]: "#4363d8",
    APRO: "#000075",
    MAD: "#469990",
    ALL: "#4363d8",
    OTHER: "#542788",
    None: None
}

MARKERS = {
    GENERAX: "o",
    RAXML: "x",
    FASTME: "d",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, APRO)]: "o",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, MAD)]: "x",
    SPEARFISH_MAP[build_ba_variant(BA_FASTME, ALL)]: "d",
    SPEARFISH_MAP[BA_FASTME]: "o",
    APRO: "o",
    MAD: "x",
    ALL: "d",
    OTHER: "o",
    None: None
}

MARKERSIZE = 2.5