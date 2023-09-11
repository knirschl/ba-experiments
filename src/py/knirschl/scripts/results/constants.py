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
