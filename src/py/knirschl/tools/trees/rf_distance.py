import os
import sys
import subprocess
import shutil
import time
from ete3 import Tree
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import paths
import utils
import fam
import metrics

def ete3_rf(tree1, tree2, unrooted = True):
  if (len(tree2.children) == 3):
    tree2.set_outgroup(tree2.children[0])
  if (len(tree1.children) == 3):
    tree1.set_outgroup(tree1.children[0])
  return tree1.robinson_foulds(tree2, unrooted_trees=unrooted, skip_large_polytomies = True, correct_by_polytomy_size = True)

def generate_scheduler_commands_file(datadir, gene_method, subst_model, output_dir):
  results_dir = os.path.join(output_dir, "results")
  os.makedirs(results_dir)
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      family_dir = fam.get_family_path(datadir, family)
      ali = fam.get_alignment(datadir, family)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      command.append("--rf")
      command.append("--tree")
      command.append(fam.build_gene_tree_path(datadir, subst_model, family, gene_method))
      command.append("--prefix")
      command.append(os.path.join(results_dir, family))
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def extract(datadir, gene_method, subst_model, output_dir):
  distances = {}  
  for family in fam.get_families_list(datadir):
    log = os.path.join(output_dir, "per_job_logs", family + "_out.txt")
    lines = open(log).readlines()
    rrf = float(lines[1].split()[-1])
    distances[family] = rrf
  return distances

def compute(datadir, gene_method, subst_model, cores):
  output_dir = fam.get_run_dir(datadir, "", "pairwise_rf_run")
  utils.reset_dir(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, gene_method, subst_model, output_dir)
  utils.run_with_scheduler(paths.raxml_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")   
  return extract(datadir, gene_method, subst_model, output_dir)


if (__name__ == "__main__"):
  if (len(sys.argv) != 5):
    print("Syntax python " + os.path.basename(__file__) + " datadir gene_method subst_model 40")
    sys.exit(1)
  datadir = sys.argv[1]
  gene_method = sys.argv[2]
  subst_model = sys.argv[3]
  cores=  sys.argv[4]
  print(compute(datadir, gene_method, subst_model, cores))