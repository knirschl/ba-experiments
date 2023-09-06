import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import paths
import utils
import fam
import metrics
import sequence_model
#import rf_cells
#import fast_rf_cells
#import extract_event_number
#import events

def get_possible_strategies():
  return ["SPR", "EVAL", "SKIP", "RECONCILE"]

def check_inputs(strategy):
  if (not (strategy in get_possible_strategies())):
    print("Unknown search strategy " + strategy)
    exit(1)


def has_multiple_sample(starting_tree):
  return "ale" in starting_tree.lower() or "multiple" in starting_tree.lower()

def get_starting_tree_path(datadir, subst_model, family, starting_tree):
  #if (has_multiple_sample(starting_tree)):
  #  return os.path.join(fam.get_family_misc_dir(datadir, family), starting_tree + "." + subst_model + "_onesample.geneTree")
  #else:
    return fam.build_gene_tree_path(datadir, subst_model, family, starting_tree)

# GeneRax does not accept tree files with multiple trees
def sample_one_starting_tree(datadir, subst_model, starting_tree):
  for family in fam.get_families_list(datadir):
    input_tree = fam.build_gene_tree_path(datadir, subst_model, family, starting_tree)
    output_tree = get_starting_tree_path(datadir, subst_model, family, starting_tree)
    tree = open(input_tree, "r").readline()
    open(output_tree, "w").write(tree)

def build_generax_families_file(datadir, starting_tree, subst_model, output):
  if (has_multiple_sample(starting_tree)):
    sample_one_starting_tree(datadir, subst_model, starting_tree)
  families_dir = os.path.join(datadir, "families")
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    print("starting gene tree " + starting_tree)
    for family in os.listdir(families_dir):
      family_path = os.path.join(families_dir, family)
      writer.write("- " + family + "\n")
      gene_tree = get_starting_tree_path(datadir, subst_model, family, starting_tree)
      if (starting_tree == "random"):
        gene_tree = "__random__"
      writer.write("starting_gene_tree = " + gene_tree + "\n")
      writer.write("alignment = " + fam.get_alignment_file(family_path) + "\n")
      writer.write("mapping = " + fam.get_mappings(datadir, family) + "\n")
      raxml_model = ""
      if (starting_tree != "random" and starting_tree != "true"):
        raxml_model = fam.get_raxml_best_model(datadir, subst_model, family)
      if (os.path.isfile(raxml_model)):
        writer.write("subst_model = " + raxml_model + "\n")
      else:
        writer.write("subst_model = " + sequence_model.get_raxml_model(subst_model) + "\n")


def build_generax_families_file_eval(datadir, subst_model, output, tree_prefix=""):
  families_dir = fam.get_families_dir(datadir)
  empty = True
  skip = 0
  with open(output, "w") as writer:
    writer.write("[FAMILIES]\n")
    for family in fam.get_families_list(datadir):
      alignment = fam.get_alignment(datadir, family)
      mapping = fam.get_mappings(datadir, family)
      raxml_model = sequence_model.get_raxml_model(subst_model)
      for tree in fam.get_gene_tree_list(datadir, family):
        if (not tree.startswith(tree_prefix)):
          continue
        scale = float(re.search(r'(\d+(?:\.\d+)?)S~G', tree)[1])
        #if (scale < 1 or scale > 4.5):
        #  continue
        treefam = family + ">" + tree.replace(".geneTree.newick", "")
        if (os.path.isfile(os.path.join(fam.get_run_dir(datadir, subst_model, "generax_eval_run"), "results", treefam, "stats.txt"))):
          # already evaluated
          skip += 1
          continue
        empty = False
        writer.write("- " + treefam + "\n")
        writer.write("starting_gene_tree = " + os.path.join(fam.get_gene_tree_dir(datadir, family), tree) + "\n")
        writer.write("alignment = " + alignment + "\n")
        writer.write("mapping = " + mapping + "\n")
        writer.write("subst_model = " + raxml_model + "\n")
    print("~~~~ Skipped", skip, "trees from getting evaluated ~~~~")
    return empty

def get_generax_command(generax_families_file, species_tree, strategy, rec_model, additional_arguments, output_dir, mode, cores):
    executable = paths.generax_exec
    old = utils.checkAndDelete("--old", additional_arguments) 
    if (mode == "gprof"):
      executable = paths.generax_gprof_exec
    elif (mode == "scalasca"):
      executable = paths.generax_scalasca_exec
    if (old):
      executable += "old"
    #generax_output = os.path.join(output_dir, "generax")
    command = []
    command.append("mpirun")
    command.append("-np")
    command.append(str(cores))
    command.append(executable)
    command.append("-f")
    command.append(generax_families_file)
    command.append("-s")
    command.append(species_tree)
    command.append("--strategy")
    command.append(strategy)
    command.append("-p")
    command.append(output_dir)
    command.extend(additional_arguments)
    return " ".join(command)

def run_generax(datadir,  subst_model,  strategy, rec_model, species_tree, generax_families_file, mode, cores, resultsdir, additional_arguments = ""):
  command = get_generax_command(generax_families_file, species_tree, strategy, rec_model, additional_arguments, resultsdir, mode, cores)
  #print(command)
  subprocess.check_call(command.split(" "), stdout = sys.stdout)


def get_mode_from_additional_arguments(additional_arguments):
  mode = "normal"
  if ("--scalasca" in additional_arguments):
    mode = "scalasca"
    additional_arguments.remove("--scalasca")
  elif ("--gprof" in additional_arguments):
    mode = "gprof"
    additional_arguments.remove("--gprof")
  return mode


def extract_trees(datadir, results_family_dir, run_name, subst_model):
  results_dir = os.path.join(results_family_dir, "reconciliations")
  for family in fam.get_families_list(datadir):
    #source = os.path.join(results_dir, family, "geneTree.newick")
    source = os.path.join(results_dir, family + "_events.newick")
    dest = fam.build_gene_tree_path_from_run(datadir, family, run_name)
    try:
      shutil.copy(source, dest)
    except:
      pass

def get_run_name(species_tree, gene_trees, subst_model, strategy, additional_arguments):
    rec_model = utils.getArg("--rec-model", additional_arguments, "UndatedDTL") 
    radius = utils.getArg("--max-spr-radius", additional_arguments, "5") 
    per_fam_rates = "--per-family-rates" in additional_arguments
    per_species_rates = "--per-species-rates" in additional_arguments
    run_name = "generax-" + strategy + "-" + rec_model + "-r" + radius
    if (per_fam_rates):
        run_name += "-famrates"
    if (per_species_rates):
        run_name += "-speciesrates"
    tc = utils.getArg("--transfer-constraint", additional_arguments, "NONE")
    if (tc == "PARENTS"):
      run_name += "-tcparent"
    if (tc == "SOFTDATED"):
      run_name += "-tcsoft"
    seed = utils.getArg("--seed", additional_arguments, None)
    if (seed != None):
      run_name += "-seed" + str(seed)
    run_name += "." + subst_model
    return run_name

def eval(results_dir, family):
  lines = open(os.path.join(results_dir, family, "stats.txt")).readlines()
  stats = {}
  logls = lines[0].split()
  stats["raxmlLogL"] = (float)(logls[0])
  stats["otherLogL"] = (float)(logls[1])
  stats["sumLogL"] = (float)(logls[0]) + (float)(logls[1])
  dlt = lines[1].split()
  stats["dup"] = (float)(dlt[3])
  stats["loss"] = (float)(dlt[4])
  #stats["transfer"] = (float)(dlt[5])
  return stats

def eval_sumLogL(results_dir, family):
  logls = open(os.path.join(results_dir, family, "stats.txt")).readlines()[0].split()
  return (float)(logls[0]) + (float)(logls[1])

def eval_and_pick(datadir, results_dir):
  best_tree = {}
  best_logL = {}
  for family in os.listdir(results_dir):
    true_family, tree = family.split(">")
    idx = 0 if "a." in tree else (1 if "m." in tree else 2)
    if (true_family not in best_tree):
      best_tree[true_family] = 3 * [""]
      best_logL[true_family] = 3 * [float("-inf")]
    logL = eval_sumLogL(results_dir, family)
    if (logL > best_logL[true_family][idx]):
      best_tree[true_family][idx] = tree
      best_logL[true_family][idx] = logL
    elif (logL == best_logL[true_family][idx]):
      # TODO if same: which one?
      continue
  with open(os.path.join(fam.get_metrics_dir(datadir), "generax_picks.txt"), "w") as writer:
    for family in best_tree:
      # pick = best_tree[family]
      for bt in best_tree[family]:
        writer.write(family + "  " + bt + "\n")
      # old_name = [f for f in fam.get_gene_tree_list(datadir, family) if f.startswith(pick)][0]
      # if (not ".generax_pick" in old_name):
      # os.rename(os.path.join(fam.get_gene_tree_dir(datadir, family), old_name), os.path.join(fam.get_gene_tree_dir(datadir, family), pick + ".generax_pick.geneTree.newick"))

def extract_events(datadir, results_family_dir, additional_arguments):
  #try:
    rec_model = utils.getArg("--rec-model", additional_arguments, "UndatedDTL")
    radius = int(utils.getArg("--max-spr-radius", additional_arguments, "5"))
    event_counts = extract_event_number.extract(results_family_dir)
    events.update_event_counts(datadir, rec_model, radius, event_counts)

def run(datadir, subst_model, strategy, species_tree, starting_tree, cores, additional_arguments, resultsdir, do_extract=True):
  run_name = utils.getAndDelete("--run", additional_arguments, None) 
  if (None == run_name):
      run_name = get_run_name(species_tree, starting_tree, subst_model, strategy, additional_arguments)
  print("Run name " + run_name)
  if (strategy == "EVAL"):
    os.makedirs(resultsdir, exist_ok=True)
  else:
    shutil.rmtree(resultsdir, True)
    os.makedirs(resultsdir)
  sys.stdout.flush()
  mode = get_mode_from_additional_arguments(additional_arguments)
  rec_model = utils.getArg("--rec-model", additional_arguments, "UndatedDTL")
  generax_families_file = os.path.join(resultsdir, "families-generax.txt")
  empty = False
  if (strategy == "EVAL"):
    empty = build_generax_families_file_eval(datadir, subst_model, generax_families_file, tree_prefix=starting_tree)
  else:
    build_generax_families_file(datadir, starting_tree, subst_model, generax_families_file)
  start = time.time()
  if (not empty):
    run_generax(datadir, subst_model, strategy, rec_model, species_tree, generax_families_file, mode, cores, resultsdir, additional_arguments)
    metrics.save_metrics(datadir, run_name, (time.time() - start), "runtimes") 
    metrics.save_metrics(datadir, run_name, (time.time() - start), "seqtimes")
  else:
    print("Skipping generax as no families are provided")
  if (strategy == "EVAL"):
    eval_and_pick(datadir, os.path.join(resultsdir, "results"))
  if (do_extract):
    extract_trees(datadir, resultsdir, run_name, subst_model)
  #print("Output in " + resultsdir)
  return resultsdir

def launch(datadir, subst_model, strategy, species_tree, starting_tree, cluster, cores, additional_arguments):
  command = ["python3"]
  command.extend(sys.argv)
  command.append("--exprun")
  dataset = os.path.basename(datadir)
  resultsdir = os.path.join("GeneRax", dataset, strategy + "_" + species_tree +  "_start_" + starting_tree, "run")
  resultsdir = utils.create_result_dir(resultsdir, additional_arguments)
  submit_path = os.path.join(resultsdir, "submit.sh")
  command.append(resultsdir)
  print(" ".join(command))
  utils.submit(submit_path, " ".join(command), cores, cluster) 
  

if (__name__ == "__main__"): 
  print("launch_generax " + str(sys.argv))
  is_run = ("--exprun" in sys.argv)
  resultsdir = ""
  if (is_run):
    resultsdir = sys.argv[-1]
    sys.argv = sys.argv[:-2]
    
  min_args_number = 8
  if (len(sys.argv) < min_args_number):
    print("Syntax error: python " + os.path.basename(__file__) + "  dataset species_tree gene_trees subst_model strategy cluster cores [additional paremeters].\n  ")
    sys.exit(1)

  datadir = os.path.normpath(sys.argv[1])
  species_tree = sys.argv[2]
  starting_tree = sys.argv[3]
  subst_model = sys.argv[4]
  strategy = sys.argv[5]
  cluster = sys.argv[6]
  cores = int(sys.argv[7])
  additional_arguments = sys.argv[min_args_number:]
  check_inputs(strategy)

  if (starting_tree == "raxml"):
    print("use raxml-ng instead of raxml please")
    exit(1)

  if (is_run):
    run(datadir, subst_model, strategy, species_tree, starting_tree, cores, additional_arguments, resultsdir)
  else:
    launch(datadir, subst_model, strategy, species_tree, starting_tree, cluster, cores, additional_arguments)
