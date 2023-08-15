import os
import re
import sys
sys.path.insert(0, 'scripts/programs')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
import launch_generax
import fam
import metrics
import rf_distance


def global_compare(root_output, replicates, tag):
    # AVERAGE OVER ALL REPLICATES
    abs_name = "rf_distance_avg-abs"
    rel_name = "rf_distance_avg-rel"
    rt_name = "runtimes"
    abs_avgs_dico = metrics.get_metrics(replicates[0], abs_name)
    rel_avgs_dico = metrics.get_metrics(replicates[0], rel_name)
    rt_avgs_dico = metrics.get_metrics(replicates[0], rt_name)
    rep_counter = 1
    best_tree_avg = 0
    best_tree_counter = 0
    for rep in replicates[1:]:
        # compute global averages
        cur_abs = metrics.get_metrics(rep, abs_name)
        cur_rel = metrics.get_metrics(rep, rel_name)
        cur_rt = metrics.get_metrics(rep, rt_name)
        for x in set(abs_avgs_dico).union(cur_abs):
            if not x in abs_avgs_dico:
                abs_avgs_dico[x] = 0
                rel_avgs_dico[x] = 0
            if not x in cur_abs:
                cur_abs[x] = 0
                cur_rel[x] = 0
            abs_avgs_dico[x] = (float(abs_avgs_dico[x]) * rep_counter + float(cur_abs[x])) / (
                        rep_counter + 1)
            rel_avgs_dico[x] = (float(rel_avgs_dico[x]) * rep_counter + float(cur_rel[x])) / (
                        rep_counter + 1)
        for x in set(rt_avgs_dico).union(cur_rt):
            if not x in rt_avgs_dico:
                rt_avgs_dico[x] = 0
            if not x in cur_rt:
                cur_rt[x] = 0
            rt_avgs_dico[x] = (float(rt_avgs_dico[x]) * rep_counter + float(cur_rt[x])) / (
                        rep_counter + 1)
        rep_counter += 1
        # get best tree distance per family and compute average over this
        for family in fam.get_families_list(rep):
            with open(os.path.join(fam.get_family_path(rep, family), "metrics",
                                   "rf_distance-rel.txt"), 'r') as reader:
                line = reader.readline()
                if (line == ''):
                    continue
                best_dist = float(line.split(" : ")[1].replace("\n", ''))
                best_tree_avg = (best_tree_avg * best_tree_counter + best_dist) / (
                            best_tree_counter + 1)
                best_tree_counter += 1
    rt_avgs_dico = {k: v for k, v in rt_avgs_dico.items() if not 'pipeline' in k}
    metrics.save_dico(root_output, abs_avgs_dico, tag + "_global__" + abs_name)
    metrics.save_dico(root_output, rel_avgs_dico, tag + "_global__" + rel_name)
    metrics.save_dico(root_output, rt_avgs_dico, tag + "_global__runtimes_avg")
    metrics.update_dico(root_output, {"best_tree_avg": best_tree_avg}, "misc")
    print("Avg best distance:", best_tree_avg)
    with open(os.path.join(fam.get_metrics_dir(root_output), "DL_global__rf_distance_avg-rel.txt"),
              'r') as reader:
        line = reader.readline()
        if (line != ''):
            split = line.split(" : ")
            return (split[0], split[1])

def collect_generax_picks(root_output, replicates, tag, compare):
    try:
        os.makedirs(os.path.join(root_output, "metrics"))
    except:
        pass
    glob = {}
    pos = 1
    try:
        with open(os.path.join(root_output, "metrics",
                               tag + "_global__rf_distance_avg-rel.txt")) as reader:
            for line in reader.readlines():
                glob[line.split(" : ")[0].lower()] = pos
                pos += 1
    except:
        glob = None
    poss = [0 for _ in range(len(glob))]
    dists = [0 for _ in range(101)]
    pick_avg_rel_dist = 1.0
    pick_ctr = 0
    with open(os.path.join(root_output, "metrics", "generax_picks.txt"), "w") as writer:
        for rep in replicates:
            seed = re.search(r'.*?seed([0-9]+)', rep)[1]
            writer.write(seed + '\n')
            with open(os.path.join(fam.get_metrics_dir(rep), "generax_picks.txt"), "r") as reader:
                for line in reader.readlines():
                    line = line[:-1]
                    writer.write(line)
                    split = line.split("  ")
                    family = split[0]
                    tree_pick = split[1] + ".geneTree.newick"
                    true_tree = fam.get_true_tree(rep, family)
                    if (glob or compare):
                        writer.write("  (")
                        if (glob):
                            pos = glob[tree_pick.lower()]
                            writer.write("pos=" + str(pos))
                            poss[pos - 1] += 1
                        if (compare):
                            dist = rf_distance.raxmlng_rf(
                                os.path.join(fam.get_gene_tree_dir(rep, family), tree_pick),
                                true_tree)
                            if (dist[0] == "abort"):
                                # print(seed, family, gt, '\n', dist)
                                writer.write(")\n")
                                continue
                            dist = dist[1]
                            # ", " if glob=True
                            writer.write(", " * (glob != None) + "dist=" + str(dist))
                            dists[int(dist * 100)] += 1
                            if (dist <= 1):  # threshold
                                pick_avg_rel_dist = (pick_avg_rel_dist * pick_ctr + dist) / (
                                            pick_ctr + 1)
                                pick_ctr += 1
                        writer.write(")")
                    writer.write('\n')

    # metrics.update_dico(root_output, {"rank-distribution.rk" + str(i + 1) : poss[i] for i in range(len(poss))}, "misc")
    # metrics.update_dico(root_output, {"distance-distribution.dist" + str(i + 1): dists[i] for i in range(len(dists))}, "misc")
    print("Rank distribution:", poss)
    print("Distance distribution:", dists)
    metrics.update_dico(root_output, {"pick_tree_avg": pick_avg_rel_dist}, "misc")
    print("Avg pick distance:", pick_avg_rel_dist)

def generax_likelihood_comp(root_output, replicates, best_avg_tree, generax_run_dir):
    stats_out = []
    with open(os.path.join(root_output, "metrics", "generax_picks.txt"), 'r') as reader:
        replicate = ""
        seed_num = 0
        for line in reader.readlines():
            if (re.match(r'[0-9]+\n', line)):
                # new dataset
                #seed = int(line.replace('\n', ''))
                replicate = replicates[seed_num]
                seed_num += 1
                continue
            # same dataset as before
            m = re.match(r'family_([0-9]+)\s+(\S*)\s+(.*)', line)
            family = m[1]
            tree_id = m[2]
            # dist == None if it didn't match
            dist = re.search(r'dist=(1(?:\.0)|0(?:\.[0-9]+))', m[3])
            if (dist):
                dist = dist[1]
            # stats = [raxmlLogLi, otherLogLi, sumLogLi, dup, loss]
            stats_picked = launch_generax.eval(os.path.join(replicate, generax_run_dir, "results"), "family_" + family + ">" + tree_id)
            try:
                stats_best_avg = launch_generax.eval(
                    os.path.join(replicate, generax_run_dir, "results"),
                    "family_" + family + ">" + best_avg_tree.replace("s~g", "S~G").replace(
                        "genetree.newick", ''))
            except FileNotFoundError as e:
                stats_best_avg = []
            rf_rel = metrics.get_metrics(fam.get_family_path(replicate, "family_" + family),
                                         "rf_distance-rel")
            # stats = [..., trueDist]
            try:
                stats_picked["rfDist"] = rf_rel[tree_id.lower() + ".genetree.newick"]
                stats_best_avg["rfDist"] = rf_rel[best_avg_tree]
            except:
                pass
            stats_out.append(
                "{} {} {} Stats_picked: {}\n".format(replicate, family, tree_id, stats_picked))
            stats_out.append(
                "{} {} {} Stats_best-on-avg: {}\n".format(replicate, family, best_avg_tree,
                                                          stats_best_avg))
    with open(os.path.join(root_output, "..", "logs", "likelihood_stats.txt"), 'w') as writer:
        for l in stats_out:
            writer.write(l)
