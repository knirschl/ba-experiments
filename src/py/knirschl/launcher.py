import os
import sys
sys.path.insert(0, 'scripts')
import paths
import utils

# SETTINGS
debug = False
cluster = "normal"
if ("basement" in os.getcwd()):
    cluster = "haswell" # + "d" if debug
elif ("fast" in os.getcwd()):
    cluster = "cascade"
reps = 100
parts = 4 # 10 should work in every case
run_filter = "full" # sim, full, fm-ba-pc, ba-pc, pc, comp 
enable_pip = True
enable_eval = parts == 1
compare_picks = True

# mistake?
if (debug):
    print("CAUTION: Running in DEBUG mode! Is this intended?")
if (parts != 1 and run_filter == "comp"):
    print("CAUTION: Multiple parts while only comparing!")
if (enable_pip and enable_eval):
    print("INFO: Computing and evaluating at the same time.")

# VALUES
benchmarks = {"BASE" : [-1],
            "SPECIES" : [15, 50, 75, 100],
            "SITES" : [50, 250, 500],
            "BRALEN" : [0.01, 0.1, 10.0, 100.0],
            "DUPLOS" : [0.0, 0.5, 2.0, 3.0]}
benchmarks_ext = {"SPECIES" : [],
            "FAM" : [10, 50, 250],
            "SITES" : [],
            "BRALEN" : [],
            "DUPLOS" : [],
            "TRA" : [],
            "DUPLOSTRA" : [],
            "POPS" : [50, 1000, 470000000]}

reps_per_part = int(reps / parts)
for bmark in benchmarks:
    for val in benchmarks[bmark]:
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
            command.append(str(int(compare_picks))) # compare generax picks
            command = " ".join(command)
            utils.submit(os.path.join(paths.output_root, "submit", "run_" + bmark + str(val) + "_part" + str(part) + ".sh"), command, 16, cluster)