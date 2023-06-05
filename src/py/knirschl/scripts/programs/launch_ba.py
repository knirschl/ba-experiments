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

def generate_scheduler_commands_file(datadir, subst_model, species_tree, ba_families_file, cores, output_dir):
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
            command.append("-s")
            command.append(species_tree)
            # CURRENTLY NOT USING FAMILIES FILE BUT STARTING NEW INSTANCE EVERY TIME
            #command.append("-f")
            #command.append(ba_families_file)
            command.append("-a")
            command.append(alignment_matrix)
            command.append("-p")
            command.append(ba_output_prefix)
            writer.write(" ".join(command) + "\n")
    return scheduler_commands_file

def extract_ba_trees(datadir, subst_model):
    families_dir = os.path.join(datadir, "families")
    valid = 0
    invalid = 0
    for family in os.listdir(families_dir):
        for miscfile in os.listdir(fam.get_family_misc_dir(datadir, family)):
            if not miscfile.startswith("ba_output." + subst_model):
                continue
            batree = os.path.join(fam.get_family_misc_dir(datadir, family), miscfile)
            tree = os.path.join(fam.get_gene_tree_dir(datadir, family), miscfile)
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
    
def build_ba_families_file(datadir, output):
    families_dir = os.path.join(datadir, "families")
    with open(output, "w") as writer:
        writer.write("[FAMILIES]\n")
        plop = 0
        for family in os.listdir(families_dir):
            family_path = os.path.join(families_dir, family)
            writer.write("- " + family + "\n")
            writer.write("alignment_matrix = " + fam.get_alignment_matrix_sorted_file(family_path) + "\n")
            #raxml_model = ""
            #if (starting_tree != "random" and starting_tree != "true"):
            #    raxml_model = fam.get_raxml_best_model(datadir, subst_model, family)
            #if (os.path.isfile(raxml_model)):
            #    writer.write("subst_model = " + raxml_model + "\n")
            #else:
            #    writer.write("subst_model = " + sequence_model.get_raxml_model(subst_model) + "\n")

def run_ba_on_families(datadir, subst_model, species_tree, cores):
    # output dir
    output_dir = fam.get_run_dir(datadir, subst_model, "ba_run")
    shutil.rmtree(output_dir, True)
    os.makedirs(output_dir)
    # config file -- UNUSED -- 
    ba_families_file = os.path.join(output_dir, "families-ba.txt")
    build_ba_families_file(datadir, ba_families_file)
    # run
    scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, species_tree, ba_families_file, cores, output_dir)
    start = time.time()
    utils.run_with_scheduler(paths.ba_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")
    # metrics
    metrics.save_metrics(datadir, fam.get_run_name("ba", subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(output_dir)
    metrics.save_metrics(datadir, fam.get_run_name("ba", subst_model), (time.time() - start) * lb, "seqtimes") 
    extract_ba_trees(datadir, subst_model)