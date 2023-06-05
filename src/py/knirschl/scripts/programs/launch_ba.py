import os
import shutil
import time
import sys
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import paths
import utils
import fam
import metrics

def generate_scheduler_commands_file(datadir, subst_model, species_tree, cores, output_dir):
    results_dir = os.path.join(output_dir, "results")
    scheduler_commands_file = os.path.join(output_dir, "commands.txt")
    with open(scheduler_commands_file, "w") as writer:
        for family in fam.get_families_list(datadir):
            ba_dir = fam.get_family_misc_dir(datadir, family)
            alignment_matrix = fam.get_alignment_matrix_sorted(datadir, family)
            try:
                os.mkdir(ba_dir)
            except:
                pass
            ba_output_prefix = os.path.join(ba_dir, "ba_output." + subst_model)
            command = []
            command.append(family)
            command.append("1")
            command.append("1")
            command.append("--alignmat") # -a
            command.append(alignment_matrix)
            command.append("--speciesmat") # -s
            command.append(species_tree)
            #command.append("--model")
            #command.append(subst_model)
            command.append("--prefix") # -p
            command.append(ba_output_prefix)
            writer.write(" ".join(command) + "\n")
    return scheduler_commands_file

def extract_ba_trees(datadir, subst_model):
    families_dir = os.path.join(datadir, "families")
    valid = 0
    invalid = 0
    for family in os.listdir(families_dir):
        batree = os.path.join(families_dir, family, "misc", "ba_output." + subst_model + ".newick")
        tree = fam.build_gene_tree_path(datadir, subst_model, family, "ba")
        if (os.path.isfile(batree) and os.stat(batree).st_size > 0):
            valid += 1
            shutil.copyfile(batree, tree)
        else:
            invalid += 1
            try:
                os.remove(tree)
            except:
                pass
        os.remove(batree) 
    print("Extracted " + str(valid) + " trees")
    if (invalid > 0):
        print("WARNING! " + str(invalid) + " trees were skipped")

def run_ba_on_families(datadir, subst_model, species_tree, cores):
    output_dir = fam.get_run_dir(datadir, subst_model, "ba_run")
    shutil.rmtree(output_dir, True)
    os.makedirs(output_dir)
    scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, species_tree, cores, output_dir)
    start = time.time()
    utils.run_with_scheduler(paths.ba_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")
    metrics.save_metrics(datadir, fam.get_run_name("ba", subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(output_dir)
    metrics.save_metrics(datadir, fam.get_run_name("ba", subst_model), (time.time() - start) * lb, "seqtimes") 
    extract_ba_trees(datadir, subst_model)