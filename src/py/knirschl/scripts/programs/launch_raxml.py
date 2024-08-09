import os
import shutil
import sys
import time

sys.path.insert(0, 'scripts')
sys.path.insert(0, os.path.join("tools", "families"))
sys.path.insert(0, 'tools/msa')
import paths
import utils
import fam
import metrics
import analyze_msa

def generate_scheduler_commands_file(datadir, subst_model, output_dir):
    results_dir = os.path.join(output_dir, "results")
    scheduler_commands_file = os.path.join(output_dir, "commands.txt")
    #gamma = False
    #sp = subst_model.split("+")
    #if (len(sp) > 1 and sp[1] == "G"):
        #fastme_model = sp[0]
        #gamma = True
    with open(scheduler_commands_file, "w") as writer:
        for family in fam.get_families_list(datadir):
            if (not analyze_msa.has_distinct_seqs(
                    fam.get_alignment_file(fam.get_family_path(datadir, family)))):
                # not enough distinct sequences
                continue
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
            command.append("--redo") # ignore checkpoints
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
            os.remove(raxmltree) 
        else:
            invalid += 1
            try:
                os.remove(tree)
            except:
                pass
    print("Extracted " + str(valid) + " trees")
    if (invalid > 0):
        print("WARNING! " + str(invalid) + " trees were skipped")

def run_raxmlng_on_families(datadir, subst_model, cores, mpi = True):
    output_dir = fam.get_run_dir(datadir, subst_model, "raxml_run")
    shutil.rmtree(output_dir, True)
    os.makedirs(output_dir)
    scheduler_commands_file = generate_scheduler_commands_file(datadir, subst_model, output_dir)
    start = time.time()
    if (mpi):
        utils.run_with_scheduler(paths.raxml_mpi_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")
    else:
        utils.run_with_scheduler(paths.raxml_exec, scheduler_commands_file, "fork", cores, output_dir, "logs.txt")   
    metrics.save_metrics(datadir, fam.get_run_name("raxml", subst_model), (time.time() - start), "runtimes") 
    lb = fam.get_lb_from_run(output_dir)
    metrics.save_metrics(datadir, fam.get_run_name("raxml", subst_model), (time.time() - start) * lb, "seqtimes") 
    extract_raxml_trees(datadir, subst_model)
    
