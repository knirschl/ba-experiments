import os
import re
import shutil
import sys
import time

sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, os.path.join("tools", "msa"))
import paths
import utils
import fam
import metrics
import analyze_msa
import msa_converter


def generate_scheduler_commands_file(datadir, subst_model, is_dna, algo, use_spr, only_mat, output_dir):
  results_dir = os.path.join(output_dir, "results")
  scheduler_commands_file = os.path.join(output_dir, "commands.txt")
  gamma = False
  sp = subst_model.split("+")
  fastme_model = subst_model
  if (len(sp) > 1 and sp[1] == "G"):
    fastme_model = sp[0]
    gamma = True
  if is_dna:
    fastme_model.removeprefix('F8')
  with open(scheduler_commands_file, "w") as writer:
    for family in fam.get_families_list(datadir):
      if (not analyze_msa.has_distinct_seqs(
              fam.get_alignment_file(fam.get_family_path(datadir, family)))):
        # not enough distinct sequences
        continue
      fastme_dir = fam.get_family_misc_dir(datadir, family)
      try:
        os.mkdir(fastme_dir)
      except:
        pass
      fastme_output = os.path.join(fastme_dir, "fastme." + subst_model + ".newick")
      fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
      command = []
      command.append(family)
      command.append("1")
      command.append("1")
      # input file (relaxed phylip format)
      command.append("-i")
      # msa
      phylip = fam.get_alignment_phylip(datadir, family)
      if (not os.path.isfile(phylip)):
        ali = fam.get_alignment(datadir, family)
        msa_converter.msa_convert(ali, phylip, None, "iphylip_relaxed")
      command.append(phylip)
      # dna or protein model
      if (is_dna):
       command.append("-d" + fastme_model)
      else:
        command.append("-p" + fastme_model)
      if (gamma):
        command.append("1.0")
      # output tree & matrix
      command.append("-o")
      command.append(fastme_output)
      command.append("-O")
      command.append(fastme_matrix)
      # distance algorithm (B: BME, I: BIONJ (def), N: NJ)
      command.append("-m")
      command.append(algo)
      command.append("-n")
      command.append("B")
      if (use_spr):
        # use spr moves
        command.append("--spr")
      if (only_mat):
        # no tree building
        command.append("-c") # only mat
        #command.append("-q") # triangular inequality correction
        command.append("-f") # precision (max=17)
        command.append("17")
      writer.write(" ".join(command) + "\n")
  return scheduler_commands_file


def extract_fastme_trees(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  valid = 0
  invalid = 0
  for family in os.listdir(families_dir):
    for miscfile in os.listdir(fam.get_family_misc_dir(datadir, family)):
      if (not (("fastme" in miscfile) and miscfile.endswith(".newick"))):
        continue
      fastmetree = os.path.join(fam.get_family_misc_dir(datadir, family), miscfile)
      tree = os.path.join(fam.get_gene_tree_dir(datadir, family), miscfile)
      #fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
      fastme_matrix = fastmetree.replace("geneTree.newick", "matrix.phy").replace("fastme.", "")
      if (os.path.isfile(fastmetree) and os.stat(fastmetree).st_size > 0):
        valid += 1
        shutil.copyfile(fastmetree, tree)
        os.remove(fastmetree)
      else:
        invalid += 1
        try:
          os.remove(tree)
          os.remove(fastme_matrix)
        except:
          pass
      try:
        os.remove(fastme_matrix + "_fastme_stat.txt")
        os.remove(fastmetree + "_fastme_stat.txt")
      except:
        pass
  print("Extracted " + str(valid) + " trees")
  if (invalid > 0):
    print("WARNING! " + str(invalid) + " trees were skipped")


def extract_fastme_mats(datadir, subst_model):
  families_dir = os.path.join(datadir, "families")
  valid = 0
  invalid = 0
  for family in os.listdir(families_dir):
    fastme_matrix = fam.get_fastme_distances(datadir, family, subst_model)
    matrix = fam.get_alignment_matrix(datadir, family)
    if (os.path.isfile(fastme_matrix) and os.stat(fastme_matrix).st_size > 0):
      valid += 1
      shutil.copyfile(fastme_matrix, matrix)
    else:
      invalid += 1
      try:
        os.remove(matrix)
        os.remove(fastme_matrix)
      except:
        pass
    try:
      os.remove(fastme_matrix)
      os.remove(fastme_matrix + "_fastme_stat.txt")
    except:
      pass
  print("Extracted " + str(valid) + " matrices")
  if (invalid > 0):
    print("WARNING! " + str(invalid) + " matrices were skipped")


def run_fastme_on_families(datadir, subst_model, is_dna, algo, use_spr, only_mat, cores):
  fastme_name = ["fastme", "fastme_mat"][only_mat]
  output_dir = fam.get_run_dir(datadir, subst_model, fastme_name + "_run")
  shutil.rmtree(output_dir, True)
  os.makedirs(output_dir)
  scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, is_dna, algo, use_spr, only_mat, output_dir)
  start = time.time()
  utils.run_with_scheduler(paths.fastme_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")
  metrics.save_metrics(datadir, fam.get_run_name(fastme_name, subst_model), (time.time() - start), "runtimes") 
  lb = fam.get_lb_from_run(output_dir)
  metrics.save_metrics(datadir, fam.get_run_name(fastme_name, subst_model), (time.time() - start) * lb, "seqtimes")
  utils.printFlush("Finished FastME, now extracting")
  if not only_mat:
    extract_fastme_trees(datadir, subst_model)
  else:
    extract_fastme_mats(datadir, subst_model)

if (__name__== "__main__"):
  max_args_number = 5
  if len(sys.argv) < max_args_number:
    print("Syntax error: python " + os.path.basename(__file__) + "  datadir subst_model is_dna cores.")
    print("Cluster can be either normal, haswell or magny")
    sys.exit(0)
  datadir = sys.argv[1]
  subst_model = sys.argv[2]
  is_dna = sys.argv[3] != "0"
  cores = int(sys.argv[4])
  run_fastme_on_families(datadir, subst_model, is_dna, False, cores)
