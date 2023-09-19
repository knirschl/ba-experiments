def build_ba_variant(method, algorithm):
    return method + " (" + algorithm + ")"

THRESHOLD = 10 # "zoom"
# trees
GENERAX = "GeneRax"
RAXML = "RAxML"
FASTME = "FastME"
BA = "BA"
BA_FASTME = "BA+FM"
OTHER = "other"
# algos
APRO = "APro" 
MAD = "MAD"
ALL = "All"
ALGO_IDS = {
    'a': APRO,
    'm': MAD,
    '+': ALL
}
# subset of trees
BA_VARIANTS = [build_ba_variant(meth, a) for a in [APRO, MAD, ALL] for meth in [BA, BA_FASTME]]
# Name : BASE
BASE = "BASE-1"
DATASETS = {"SPECIES" : 25,
            "SITES" : 100,
            "BRALEN" : 1.0,
            "DUPLOS" : 1.0
            }

SPEARFISH_MAP = {build_ba_variant(m, a): a.lower() + ["FM", "NJ"][m == BA] for a in [APRO, MAD, ALL] for m in [BA, BA_FASTME]}
SPEARFISH_MAP["BA+FM"] = "Spearfish"
SPEARFISH_MAP[GENERAX] = GENERAX
SPEARFISH_MAP[RAXML] = RAXML
SPEARFISH_MAP[FASTME] = FASTME