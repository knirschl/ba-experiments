from constants import *

COLORS = {
    GENERAX: "mediumblue",
    RAXML: "deepskyblue",
    FASTME: "darkviolet",
    build_ba_variant(BA, APRO): "lightgrey",
    build_ba_variant(BA, MAD): "darkgrey",
    build_ba_variant(BA, ALL): "dimgrey",
    build_ba_variant(BA_FASTME, APRO): "limegreen",
    build_ba_variant(BA_FASTME, MAD): "darkgreen",
    build_ba_variant(BA_FASTME, ALL): "yellowgreen",
    OTHER: "brown",
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
    OTHER: "o",
    None: None
}

MARKERSIZE = 2.5