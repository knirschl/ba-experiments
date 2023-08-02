import os
import sys
from threading import Thread
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import fam
import metrics
import rf_distance


def ttask(datadir, tree, avg_abs_dico, avg_rel_dico, families_dico):
    fam_counter = 0
    for family in fam.get_families_list(datadir):
        true_tree = fam.get_true_tree(datadir, family)
        abs_path_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
        if (not os.path.isfile(abs_path_tree)):
            picked_abs_path_tree = abs_path_tree.replace(".geneTree.newick", ".generax_pick.geneTree.newick")
            if (not os.path.isfile(picked_abs_path_tree)):
                #print("Tree file does not exist:", family, tree)
                continue
            abs_path_tree = picked_abs_path_tree
            #print("Picked", tree, "  (path", abs_path_tree, ")")
        # CALCULATIONS
        dist_abs, dist_rel = rf_distance.raxmlng_rf(abs_path_tree, true_tree)
        if dist_abs == "abort":
            # distance error
            #print("I've catched an error, maybe look into this?\n", dist_rel, "\n")
            continue
        # single distance
        # This should be thread-safe
        families_dico[family][tree] = dist_rel
        # method average
        if not tree in avg_abs_dico:
            # if tree (= method) not in dico, add
            avg_abs_dico[tree] = 0
            avg_rel_dico[tree] = 0
        avg_abs_dico[tree] = (avg_abs_dico[tree] * fam_counter + dist_abs) / (fam_counter + 1)
        avg_rel_dico[tree] = (avg_rel_dico[tree] * fam_counter + dist_rel) / (fam_counter + 1)
        fam_counter += 1

def compare_all(datadir):
    avg_abs_dico = {}
    avg_rel_dico = {}
    families_dico = {}
    for family in fam.get_families_list(datadir):
        families_dico[family] = {}
    threads = []
    for tree in fam.get_gene_tree_list(datadir, fam.get_families_list(datadir)[0]):
        if (tree == fam.get_true_gene_tree_name()):
            continue
        tree = tree.replace(".generax_pick", "")
        fam_thread = Thread(target=ttask, args=(datadir, tree, avg_abs_dico, avg_rel_dico, families_dico))
        threads.append(fam_thread)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # save method average distance
    metrics.save_dico(datadir, avg_abs_dico, "rf_distance_avg-abs")
    metrics.save_dico(datadir, avg_rel_dico, "rf_distance_avg-rel")
    for family in fam.get_families_list(datadir):
        metrics.save_dico(fam.get_family_path(datadir, family), families_dico[family], "rf_distance-rel")