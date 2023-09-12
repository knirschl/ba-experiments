# Common variables used in the scripts

import os

def get_parent_path(path):
  return os.path.abspath(os.path.join(path, os.pardir))

def python():
  return "python"

def python3():
  return "python3"

## ---- DIRECTORIES ---- 
root = get_parent_path(get_parent_path(os.path.realpath(__file__)))

# scripts
scripts_root = os.path.join(root, "scripts")
# tools
tools_root = os.path.join(root, "tools")

# github
ba_github_root = "/hits/basement/cme/knirsch/github/"
code_github_root = os.path.join(ba_github_root, "BA-Code")

# results
output_root = os.path.join(code_github_root, "output")
families_datasets_root = os.path.join(output_root, "families")

# ---- EXTERNAL PROGRAMS ----
programs_root = ba_github_root #os.path.join(ba_github_root, "resources", "tools")
# SimPhy
simphy_root = os.path.join(programs_root, "SimPhy")
simphy_exec = os.path.join(simphy_root, "simphy")
simphy_indelible_wrapper_exec = os.path.join(simphy_root, "INDELIble_wrapper.pl")
# RAxML-NG
raxml_root = os.path.join(programs_root, "RAxML-NG")
raxml_exec = os.path.join(raxml_root, "bin", "raxml-ng")
raxml_mpi_exec = os.path.join(raxml_root, "bin", "raxml-ng-mpi")
# GeneRax
generax_exec = os.path.join(programs_root, "GeneRax", "build", "bin", "generax")
# FastME
fastme_exec = os.path.join(programs_root, "FastME", "bin",  "fastme-2.1.6.2-linux64")
# Bachelor thesis
ba_exec = os.path.join(code_github_root, "build", "BA")
# MPI
mpischeduler_exec = os.path.join(programs_root, "MPIScheduler", "bin", "mpi-scheduler")
# constants
mpi_scheduler_heuristic = "--split-scheduler"
historic = os.path.join(root, "historic.txt")
