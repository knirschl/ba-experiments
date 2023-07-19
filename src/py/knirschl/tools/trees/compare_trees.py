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

def merge_dicts(dicts): # https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries (havoc_method)
    def reducer(accumulator, element):
        for key, value in element.items():
            accumulator[key] = accumulator.get(key, 0) + value
        return accumulator
    return reduce(reducer, dicts, {})

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

def ttask(datadir, family, avg_abs_dicos, avg_rel_dicos, fam_counter):
    true_tree = fam.get_true_tree(datadir, family)
    for tree in os.listdir(fam.get_gene_tree_dir(datadir, family)):
        abs_path_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
        if (abs_path_tree == true_tree):
            continue
        # CALCULATIONS
        dist_abs, dist_rel = rf_compare(abs_path_tree, true_tree)
        if dist_abs == "abort":
            # distance error
            #print("I've catched an error, maybe look into this?\n", dist_rel, "\n")
            continue
        # method average
        if not tree in avg_abs_dicos[fam_counter]:
            # if tree (= method) not in dico, add
            avg_abs_dicos[fam_counter][tree] = 0
            avg_rel_dicos[fam_counter][tree] = 0
        avg_abs_dicos[fam_counter][tree] = (avg_abs_dicos[fam_counter][tree] * fam_counter + dist_abs) / (fam_counter + 1)
        avg_rel_dicos[fam_counter][tree] = (avg_rel_dicos[fam_counter][tree] * fam_counter + dist_rel) / (fam_counter + 1)

        # SAVING TO FILE
        # save single distance
        # TODO !! RACE CONDITION WHILE WRITING !!
        #metrics.save_metrics(datadir, make_key(family, tree), dist_abs, "rf_distance-abs")
        #metrics.save_metrics(datadir, make_key(family, tree), dist_rel, "rf_distance-rel")

def compare_all(datadir):
    rf_avg_abs_dicos = [{} for _ in range(len(fam.get_families_list(datadir)))]
    rf_avg_rel_dicos = [{} for _ in range(len(fam.get_families_list(datadir)))]
    fam_counter = 0
    threads = []
    for family in fam.get_families_list(datadir):
        fam_thread = Thread(target=ttask, args=(datadir, family, rf_avg_abs_dicos, rf_avg_rel_dicos, fam_counter))
        fam_counter += 1
        threads.append(fam_thread)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    #print(rf_avg_abs_dicos)
    avg_abs_dico = merge_dicts(rf_avg_abs_dicos)
    #print(avg_abs_dico)
    avg_rel_dico = merge_dicts(rf_avg_rel_dicos)
    # save method average distance
    metrics.save_dico(datadir, avg_abs_dico, "rf_distance_avg-abs")
    metrics.save_dico(datadir, avg_rel_dico, "rf_distance_avg-rel")