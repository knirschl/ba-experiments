# Common variables used in the scripts

import os

def get_parent_path(path):
    return os.path.abspath(os.path.join(path, os.pardir))

def python():
    return "python"

def python3():
    return "python3"

## ---- DIRECTORIES ---- 
# local
balin_root = "/home/balin/Documents/Programming/HITS"
fili_root = "/home/fili/Documents/Programming/HITS"

# github
cluster_basement_root = "/hits/basement/cme/knirsch/github/"

cwd = os.getcwd()
if (cwd.startswith("/hits")):
    root = cluster_basement_root
    programs_root = root
elif (cwd.startswith("/home/balin")):
    root = balin_root
    programs_root = os.path.join(root, "Bioinformatics")
elif (cwd.startswith("/home/fili")):
    root = fili_root
    programs_root = os.path.join(root, "Bioinformatics")
else:
    root = cwd
    programs_root = root

# results
datasets_root = os.path.join(root, "datasets")
families_datasets_root = os.path.join(datasets_root, "families")

# ---- EXTERNAL PROGRAMS ----
# SimPhy
simphy_root = os.path.join(programs_root, "SimPhy")
simphy_exec = os.path.join(simphy_root, "bin", "simphy")
simphy_indelible_wrapper_exec = os.path.join(simphy_root, "bin", "INDELIble_wrapper.pl")
# RAxML-NG
raxml_root = os.path.join(programs_root, "RAxML-NG")
raxml_exec = os.path.join(raxml_root, "bin", "raxml-ng")
raxml_mpi_exec = os.path.join(raxml_root, "bin", "raxml-ng-mpi")
# GeneRax
generax_exec = os.path.join(programs_root, "GeneRax", "build", "bin", "generax")
# FastME
fastme_exec = os.path.join(programs_root, "FastME", "bin",  "fastme")
# Spearfish
spearfish_exec = os.path.join(root, "Spearfish", "build", "spearfish")
spearfish_wrapper_py = os.path.join(root, "Spearfish", "src", "py", "wrapper.py")
# MPI
mpischeduler_exec = os.path.join(programs_root, "MPIScheduler", "bin", "mpi-scheduler")
# constants
mpi_scheduler_heuristic = "fork"  # "split"
historic = os.path.join(root, "historic.txt")
