import os
import sys
import subprocess
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import paths
import fam
import metrics

def make_key(family, tree):
    return family + '\t' + tree 

def compare(tree1, tree2):
    command = []
    command.append(paths.raxml_exec)
    command.append("--rf")
    command.append(tree1 + "," + tree2)
  
    out = subprocess.check_output(command).decode("utf-8")
    lines = out.split("\n")
    rf_abs = lines[0].split(" ")[-1]
    rf_rel = lines[1].split(" ")[-1]
    print(rf_abs, rf_rel)
    return [float(rf_abs), float(rf_rel)]

def compare_all(datadir):
    for family in fam.get_families_list(datadir):
        raxml_dir = fam.get_family_misc_dir(datadir, family)
        try:
            os.mkdir(raxml_dir)
        except:
            pass
        raxml_output = os.path.join(raxml_dir, "distances.")
        true_tree = fam.get_true_tree(datadir, family)
        for tree in os.listdir(fam.get_gene_tree_dir(datadir, family)):
            abs_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
            if (abs_tree == true_tree):
                continue
            dist = compare(abs_tree, true_tree)
            metrics.save_metrics(datadir, make_key(family, tree), dist[0], "tree_distance-abs")
            metrics.save_metrics(datadir, make_key(family, tree), dist[1], "tree_distance-rel")