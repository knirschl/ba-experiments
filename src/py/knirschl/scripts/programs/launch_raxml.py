import os
import sys
import subprocess
import shutil
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
import paths
import utils
import fam
import metrics

def generate_scheduler_commands_file(datadir, subst_model, output_dir):
    results_dir = os.path.join(output_dir, "results")
    scheduler_commands_file = os.path.join(output_dir, "commands.txt")
    #gamma = False
    #sp = subst_model.split("+")
    #fastme_model = subst_model
    #if (len(sp) > 1 and sp[1] == "G"):
        #fastme_model = sp[0]
        #gamma = True
    with open(scheduler_commands_file, "w") as writer:
        for family in fam.get_families_list(datadir):
            raxml_dir = fam.get_family_misc_dir(datadir, family)
            alignment = fam.get_alignment(datadir, family)
            try:
                os.mkdir(raxml_dir)
            except:
                pass
            raxml_output_prefix = os.path.join(raxml_dir, "raxml_output." + subst_model)
            command = []
            command.append(family)
            command.append("1")
            command.append("1")
            command.append("--msa")
            command.append(alignment)
            command.append("--model")
            command.append(subst_model)
            command.append("--prefix")
            command.append(raxml_output_prefix)
            writer.write(" ".join(command) + "\n")
    return scheduler_commands_file

def extract_raxml_trees(datadir, subst_model):
    #families_dir = fam.get_families_dir()
    families_dir = os.path.join(datadir, "families")
    valid = 0
    invalid = 0
    for family in fam.get_families_list(datadir):
        raxmltree = os.path.join(families_dir, family, "misc", "raxml_output." + subst_model + ".raxml.bestTree")
        tree = fam.build_gene_tree_path(datadir, subst_model, family, "raxml")
        if (os.path.isfile(raxmltree) and os.stat(raxmltree).st_size > 0):
            valid += 1
            shutil.copyfile(raxmltree, tree)
        else:
            invalid += 1
            try:
                os.remove(tree)
            except:
                pass
        os.remove(raxmltree) 
    print("Extracted " + str(valid) + " trees")
    if (invalid > 0):
        print("WARNING! " + str(invalid) + " trees were skipped")

def run_raxmlng_on_families(datadir, subst_model, cores, mpi = False):
    if (mpi):
        print("The MPI version of RAxML-NG is currently not supported!")
        return
    output_dir = fam.get_run_dir(datadir, subst_model, "raxml_run")
    shutil.rmtree(output_dir, True)
    os.makedirs(output_dir)
    scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, output_dir)
    start = time.time()
    utils.run_with_scheduler(paths.raxml_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")   
    metrics.save_metrics(datadir, fam.get_run_name("raxml", subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(output_dir)
    metrics.save_metrics(datadir, fam.get_run_name("raxml", subst_model), (time.time() - start) * lb, "seqtimes") 
    extract_raxml_trees(datadir, subst_model)
    