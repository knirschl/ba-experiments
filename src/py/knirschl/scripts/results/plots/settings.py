from constants import *

COLORS = {
    GENERAX: "#e6194B",
    RAXML: "#800000",
    FASTME: "#f58231",
    build_ba_variant(BA, APRO): "lightgrey",
    build_ba_variant(BA, MAD): "darkgrey",
    build_ba_variant(BA, ALL): "dimgrey",
    build_ba_variant(BA_FASTME, APRO): "#000075",
    build_ba_variant(BA_FASTME, MAD): "#469990",
    build_ba_variant(BA_FASTME, ALL): "#4363d8",
    BA_FASTME: "#000075",
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
    build_ba_variant(BA, APRO): "o",
    build_ba_variant(BA, MAD): "x",
    build_ba_variant(BA, ALL): "d",
    build_ba_variant(BA_FASTME, APRO): "o",
    build_ba_variant(BA_FASTME, MAD): "x",
    build_ba_variant(BA_FASTME, ALL): "d",
    BA_FASTME: "o",
    APRO: "o",
    MAD: "x",
    ALL: "d",
    OTHER: "o",
    None: None
}

MARKERSIZE = 2.5