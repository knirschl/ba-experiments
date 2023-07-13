import os
import subprocess
import sys
from threading import Thread

sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import paths
import fam
import metrics

def make_key(family, tree):
    return family + '\t' + tree 

def rf_compare(tree1, tree2):
    command = []
    command.append(paths.raxml_exec)
    command.append("--rf")
    command.append(tree1 + "," + tree2)
    try:
        out = subprocess.check_output(command).decode("utf-8")
    except subprocess.CalledProcessError as error:
        return ["abort", error]
    lines = out.split("\n")
    rf_abs = lines[0].split(" ")[-1]
    rf_rel = lines[1].split(" ")[-1]
    print(rf_abs, rf_rel)
    return [float(rf_abs), float(rf_rel)]

def ttask(datadir, family, method_avg_abs_dico, method_avg_rel_dico, fam_counter):
    true_tree = fam.get_true_tree(datadir, family)
    for tree in os.listdir(fam.get_gene_tree_dir(datadir, family)):
        abs_path_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
        if (abs_path_tree == true_tree):
            continue
        # CALCULATIONS
        dist_abs, dist_rel = rf_compare(abs_path_tree, true_tree)
        if dist_abs == "abort":
            # distance error
            print("I've catched an error, maybe look into this?\n", dist_rel, "\n")
            continue
        # method average
        if not tree in method_avg_abs_dico:
            # if tree (= method) not in dico, add
            method_avg_abs_dico[tree] = 0
            method_avg_rel_dico[tree] = 0
        # TODO !! RACE CONDITION WHILE SETTING !!
        method_avg_abs_dico[tree] = (method_avg_abs_dico[tree] * fam_counter + dist_abs) / (fam_counter + 1)
        method_avg_rel_dico[tree] = (method_avg_rel_dico[tree] * fam_counter + dist_rel) / (fam_counter + 1)

        # SAVING TO FILE
        # save single distance
        # TODO !! RACE CONDITION WHILE WRITING !!
        metrics.save_metrics(datadir, make_key(family, tree), dist_abs, "rf_distance-abs")
        metrics.save_metrics(datadir, make_key(family, tree), dist_rel, "rf_distance-rel")

def compare_all(datadir):
    method_avg_abs_dico = {}
    method_avg_rel_dico = {}
    fam_counter = 0
    threads = []
    for family in fam.get_families_list(datadir):
        fam_thread = Thread(target=ttask, args=(datadir, family, method_avg_abs_dico, method_avg_rel_dico, fam_counter))
        fam_counter += 1
        threads.append(fam_thread)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # save method average distance
    metrics.save_dico(datadir, method_avg_abs_dico, "rf_distance_avg-abs")
    metrics.save_dico(datadir, method_avg_rel_dico, "rf_distance_avg-rel")