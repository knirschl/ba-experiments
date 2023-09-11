import os
import sys
sys.path.insert(0, 'scripts')
import paths
import utils

# SETTINGS
debug = False
reps = 100
parts = 1 # 10 should work in every case
run_filter = "full"  
enable_pip = False
enable_eval = run_filter != 's' and parts == 1
bm_set = 0 # 0: benchmarks, 1: benchmarks_ext, 2: both

# mistake?
if (debug):
    print("CAUTION: Running in DEBUG mode! Is this intended [y,n]?")
    if input() == "n":
        debug = False
        print("Disabeled debug mode")
if (parts != 1 and run_filter == "comp"):
    print("CAUTION: Multiple parts while only comparing! Change [n,p,f]?")
    match input():
        case 'p':
            parts = 1
            print("Set parts to '1'")
        case 'f':
            run_filter = "full"
            print("Set run_filter to 'full'")
if (enable_pip and enable_eval):
    print("INFO: Computing and evaluating at the same time. Change [n,c,e]?")
    match input():
        case 'e':
            enable_pip = False
            print("Disabled pipeline")
        case 'c':
            enable_eval = False
            print("Disabled evaluation")
if (enable_eval and (parts != 1 or run_filter != "full")):
    print("INFO: Evaluating only part of benchmark. Is this intended [y,n]?")
    if input() == 'n':
        parts = 1
        run_filter = "full"
        print("Enabled evaluation of whole benchmark")

cluster = "normal"
if ("basement" in os.getcwd()):
    cluster = "haswell"
elif ("fast" in os.getcwd()):
    cluster = "cascade"
cluster + 'd' * debug

# VALUES
benchmarks = {"BASE" : [-1],
            "SPECIES" : [15, 50, 75, 100],
            "SITES" : [50, 250, 500],
            "BRALEN" : [0.01, 0.1, 10.0, 100.0],
            "DUPLOS" : [0.0, 0.5, 2.0, 3.0]}
benchmarks_ext = {"SPECIES" : [250, 500],
            "SITES" : [1000],
            "BRALEN" : [100.0],
            "DUPLOS" : [],
            "TRA" : [0.5, 1.0, 2.0, 3.0],
            "DUPLOSTRA" : [],#"0.5,0.5", "1.0,1.0", "0.5,1.0", "1.0,0.5", "2.0,2.0"],
            "POP" : [10000000, 100000000, 1000000000]} # e7, e8, e9

bm = {}
if bm_set == 0:
    bm = benchmarks
elif bm_set == 1:
    bm = benchmarks_ext
elif bm_set == 2:
    for key in set(benchmarks.keys()).union(set(benchmarks_ext.keys())):
        bm[key] = benchmarks[key] if key in benchmarks else []
        if (key in benchmarks_ext and (vs := [v for v in benchmarks_ext[key] if v not in bm[key]]) != []):
            bm[key].append(vs)
reps_per_part = int(reps / parts)

for bmark in bm:
    for val in bm[bmark]:
        for part in range(parts):
            command = []
            if (cluster == "normal"):
                command.append(os.path.join(paths.code_github_root, ".venv", "bin", "python"))
            else:
                command.append(paths.python())
            command.append(os.path.join(paths.scripts_root, "pipeline.py"))
            command.append(str(part * reps_per_part)) # start replicate of part
            command.append(str(reps_per_part)) # replicates per part  (all if only 1 part)
            command.append(bmark) # tag / dataset
            command.append(str(val)) # value to change
            command.append(run_filter) # run filter (sim, full, fm-ba-pc, ba-pc, pc, comp)
            command.append(str(int(enable_pip))) # enable pipeline
            command.append(str(int(enable_eval))) # enable evaluation
            command.append("1") # compare generax pick distances
            command = " ".join(command)
            utils.submit(os.path.join(paths.output_root, "submit", "run_" + bmark + str(val) + "_part" + str(part) + ".sh"), command, 16, cluster)
