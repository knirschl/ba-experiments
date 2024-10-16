import os
import sys
sys.path.insert(0, 'scripts')
import paths
import utils

# SETTINGS
debug = False # Give debug priority to all runs
reps = 100 # Number of datasets to work on
parts = 20 # Divide each dataset in `parts` parts ## 10 should work in every case
run_filter = "c_PLG" # s,r,g,f,b,c
bm_sim_set = 0 # 0: benchmarks_sim, 1: benchmarks_sim_ext, 2: both
datasets = "real" # sim: simulated datasets, real: real datasets
evaluate = False # False: infer trees, True: evaluate tree distances

# DATASETS
# simulated
benchmarks_sim = {"BASE" : [-1],
            "SPECIES" : [15, 50, 75, 100],
            "SITES" : [50, 250, 500],
            "BRALEN" : [0.01, 0.1, 10.0, 100.0],
            "DUPLOS" : [0.0, 0.5, 2.0, 3.0]}
benchmarks_sim_ext = {"SPECIES" : [250, 500],
            "SITES" : [1000],
            "BRALEN" : [100.0],
            "DUPLOS" : [],
            "TRA" : [0.5, 1.0, 2.0, 3.0],
            "DUPLOSTRA" : [],#"0.5,0.5", "1.0,1.0", "0.5,1.0", "1.0,0.5", "2.0,2.0"],
            "POP" : [10000000, 100000000, 1000000000]} # e7, e8, e9

bm_sim = {}
if bm_sim_set == 0:
    bm_sim = benchmarks_sim
elif bm_sim_set == 1:
    bm_sim = benchmarks_sim_ext
elif bm_sim_set == 2:
    for key in set(benchmarks_sim.keys()).union(set(benchmarks_sim_ext.keys())):
        bm_sim[key] = benchmarks_sim[key] if key in benchmarks_sim else []
        if (key in benchmarks_sim_ext and (vs := [v for v in benchmarks_sim_ext[key] if v not in bm_sim[key]]) != []):
            bm_sim[key].append(vs)
reps_per_part = int(reps / parts)

# real
benchmarks_real = ["cyano_empirical"]
#benchmarks_real = ["generax_data/cyano_empirical"]

# mistake?
if (debug):
    print("CAUTION: Running in DEBUG mode! Is this intended [y,n]?")
    if input() == "n":
        debug = False
        print("Disabeled debug mode")
if (parts != 1 and run_filter == "c"):
    print("CAUTION: Multiple parts while only comparing! Change [n,p,f]?")
    match input():
        case 'p':
            parts = 1
            print("Set parts to '1'")
        case 'f':
            run_filter = "full"
            print("Set run_filter to 'full'")
if (evaluate):
    print("INFO: Everything else finished? [y, n]?")
    if input() == 'n':
        exit()

cluster = "normal"
if ("basement" in os.getcwd()):
    cluster = "haswell"
elif ("fast" in os.getcwd()):
    cluster = "cascade"
cluster += 'd' * debug

if datasets == "sim":
    for bmark in bm_sim:
        for val in bm_sim[bmark]:
            for part in range(parts):
                command = []
                if (cluster == "normal"):
                    command.append(os.path.join(paths.root, "ba-experiments", ".venv", "bin", "python"))
                else:
                    command.append(paths.python())
                command.append(os.path.join(paths.programs_root, "ba-experiments", "src", "py", "knirschl", "scripts", "pipeline.py"))
                command.append(["sim", "eval"][evaluate])
                command.append(bmark) # tag / dataset
                command.append(str(val)) # value to change
                if not evaluate:
                    command.append(str(part)) # rank / id of part
                    command.append(str(reps_per_part)) # replicates per part  (all if only 1 part)
                    command.append(run_filter) # run filter
                command = " ".join(command)
                utils.submit(os.path.join(paths.datasets_root, "submit", "run_" + bmark + str(val) + "_part" + str(part) + ".sh"), command, 16, cluster)
elif datasets == "real":
    for bmark in benchmarks_real:
        command = []
        if (cluster == "normal"):
            command.append(os.path.join(paths.root, "Spearfish", ".venv", "bin", "python"))
        else:
            command.append(paths.python())
        command.append(os.path.join(paths.root, "ba-experiments", "src", "py", "knirschl", "scripts", "pipeline.py"))
        command.append(["real", "eval"][evaluate])
        command.append(bmark)
        if not evaluate:
            command.append("0") # unused rank
            command.append("0") # unused n_families
            command.append(run_filter) # run filter
        command = " ".join(command)
        utils.submit(os.path.join(paths.datasets_root, "submit", "run_" + bmark + ".sh"), command, 16, cluster)

