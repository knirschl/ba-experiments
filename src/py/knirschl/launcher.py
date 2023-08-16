import os
import sys
sys.path.insert(0, 'scripts')
import paths
import utils

# SETTINGS
cluster = "normal"
if ("basement" in os.getcwd()):
    cluster = "haswell" # + "d" if debug
elif ("fast" in os.getcwd()):
    cluster = "cascade"
reps = 1
run_filter = "full"
enable_pip = True
compare_picks = True
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
            "POP" : [50, 1000, 470000000]}

for bmark in benchmarks:
    for val in benchmarks[bmark]:
        command = []
        if (cluster == "normal"):
            command.append(os.path.join(paths.code_github_root, ".venv", "bin", "python"))
        else:
            command.append(paths.python())
        command.append(os.path.join(paths.scripts_root, "pipeline.py"))
        command.append(str(reps)) # replicates
        command.append(bmark) # tag / dataset
        command.append(str(val)) # value to change
        command.append(run_filter) # run filter (sim, full, fm-ba-pc, ba-pc, pc, comp)
        command.append(str(int(enable_pip))) # enable pipeline
        command.append(str(int(compare_picks))) # compare generax picks
        command = " ".join(command)
        utils.submit(os.path.join(paths.output_root, "submit", "run_" + bmark + str(val) + ".sh"), command, 16, cluster)
        exit() # TODO REMOVE
