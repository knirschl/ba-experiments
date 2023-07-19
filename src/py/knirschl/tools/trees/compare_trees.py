import os
import subprocess
import sys
from threading import Thread
from functools import reduce
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import paths
import fam
import metrics

"""
TODO DEPRECATED
def merge_dicts(dicts): # https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries (havoc_method)
    def reducer(accumulator, element):
        for key, value in element.items():
            accumulator[key] = accumulator.get(key, 0) + value
        return accumulator
    return reduce(reducer, dicts, {})
"""

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
    #print(rf_abs, rf_rel)
    return [float(rf_abs), float(rf_rel)]

def ttask(datadir, tree, avg_abs_dico, avg_rel_dico):
    fam_counter = 0
    for family in fam.get_families_list(datadir):
        true_tree = fam.get_true_tree(datadir, family)
        abs_path_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
        # CALCULATIONS
        dist_abs, dist_rel = rf_compare(abs_path_tree, true_tree)
        if dist_abs == "abort":
            # distance error
            #print("I've catched an error, maybe look into this?\n", dist_rel, "\n")
            continue
        # method average
        if not tree in avg_abs_dico:
            # if tree (= method) not in dico, add
            avg_abs_dico[tree] = 0
            avg_rel_dico[tree] = 0
        avg_abs_dico[tree] = (avg_abs_dico[tree] * fam_counter + dist_abs) / (fam_counter + 1)
        avg_rel_dico[tree] = (avg_rel_dico[tree] * fam_counter + dist_rel) / (fam_counter + 1)
        fam_counter += 1
        
        # SAVING TO FILE
        # save single distance
        # TODO !! RACE CONDITION WHILE WRITING !!
        #metrics.save_metrics(datadir, make_key(family, tree), dist_abs, "rf_distance-abs")
        #metrics.save_metrics(datadir, make_key(family, tree), dist_rel, "rf_distance-rel")

def compare_all(datadir):
    avg_abs_dico = {}
    avg_rel_dico = {}
    threads = []
    for tree in fam.get_gene_tree_list(datadir, fam.get_families_list(datadir)[0]):
        if (tree == fam.get_true_gene_tree_name()):
            continue
        fam_thread = Thread(target=ttask, args=(datadir, tree, avg_abs_dico, avg_rel_dico))
        threads.append(fam_thread)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # save method average distance
    metrics.save_dico(datadir, avg_abs_dico, "rf_distance_avg-abs")
    metrics.save_dico(datadir, avg_rel_dico, "rf_distance_avg-rel")