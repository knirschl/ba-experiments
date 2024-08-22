import os
import re
import glob
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
    abs_avgs_dico = {}
    rel_avgs_dico = {}
    rt_avgs_dico = {}
    best_tree_avg = 0
    best_tree_counter = 0
    for rep in replicates:
        # compute global averages
        cur_abs = metrics.get_metrics(rep, abs_name)
        cur_rel = metrics.get_metrics(rep, rel_name)
        cur_rt = metrics.get_metrics(rep, rt_name)
        for x in cur_abs:
            if not x in abs_avgs_dico:
                abs_avgs_dico[x] = []
                rel_avgs_dico[x] = []
            abs_avgs_dico[x].append(float(cur_abs[x]))
            rel_avgs_dico[x].append(float(cur_rel[x]))
        for x in cur_rt:
            if not x in rt_avgs_dico:
                rt_avgs_dico[x] = []
            rt_avgs_dico[x].append(float(cur_rt[x]))
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
    # compute averages
    for x in abs_avgs_dico:
        abs_avgs_dico[x] = sum(abs_avgs_dico[x]) / len(abs_avgs_dico[x])
        rel_avgs_dico[x] = sum(rel_avgs_dico[x]) / len(rel_avgs_dico[x])
    for x in rt_avgs_dico:
        if (len(rt_avgs_dico[x]) > 2):
            rt_avgs_dico[x].remove(max(rt_avgs_dico[x]))
            rt_avgs_dico[x].remove(min(rt_avgs_dico[x]))
        rt_avgs_dico[x] = sum(rt_avgs_dico[x]) / len(rt_avgs_dico[x])
    # write out
    rt_avgs_dico = {k: v for k, v in rt_avgs_dico.items() if not 'pipeline' in k}
    metrics.save_dico(root_output, abs_avgs_dico, tag + "_global__" + abs_name)
    metrics.save_dico(root_output, rel_avgs_dico, tag + "_global__" + rel_name)
    metrics.save_dico(root_output, rt_avgs_dico, tag + "_global__runtimes_avg")
    metrics.update_dico(root_output, {"best_tree_avg": best_tree_avg}, "misc")
    print("Avg best distance:", best_tree_avg)
    with open(os.path.join(fam.get_metrics_dir(root_output), tag + "_global__rf_distance_avg-rel.txt"),
              'r') as reader:
        line = reader.readline()
        if (line != ''):
            split = line.split(" : ")
            return (split[0], split[1])

def collect_generax_picks(root_output, replicates, tag, compare):
    glob_avg_rel_name = "_global__rf_distance_avg-rel.txt"
    try:
        os.makedirs(os.path.join(root_output, "metrics"))
    except:
        pass
    glob_pos = {}
    try:
        for f in glob.glob(fam.get_metrics_dir(root_output) + "/" + tag + "*" + glob_avg_rel_name):
            pos = 1
            with open(f, 'r') as reader:
                for line in reader.readlines():
                    tree = line.split(" : ")[0].lower()
                    if not tree in glob_pos:
                        glob_pos[tree] = []
                    glob_pos[tree].append(pos) 
                    pos += 1
        for key in glob_pos:
            glob_pos[key] = sum(glob_pos[key]) / len(glob_pos[key])
    except:
        glob_pos = None
    poss = [0 for _ in range(len(glob_pos))]
    dists = [3* [0] for _ in range(101)]
    pick_avg_rel_dist = 3 * [1.0]
    pick_ctr = 3 * [0]
    with open(os.path.join(root_output, "metrics", tag + "_global__generax_picks.txt"), "w") as writer:
        for rep in replicates:
            #seed = re.search(r'.*?seed([0-9]+)', rep)[1]
            #writer.write(seed + '\n')
            with open(os.path.join(fam.get_metrics_dir(rep), "generax_picks.txt"), "r") as reader:
                for line in reader.readlines():
                    line = line[:-1]
                    split = line.split("  ")
                    family = split[0]
                    tree_pick = split[1] + ".geneTree.newick"
                    if (tree_pick == ".geneTree.newick" or tree_pick == ''):
                        print("Missing tree in", rep, family)
                        continue
                    writer.write(line)
                    cleaned = split[1].replace("ba.", "").replace("fastme.", "")
                    idx = 0 if "a." in cleaned else (1 if "m." in cleaned else 2)
                    true_tree = fam.get_true_tree(rep, family)
                    if (glob_pos or compare):
                        writer.write("  (")
                        if (glob_pos):
                            pos = glob_pos[tree_pick.lower()]
                            writer.write("pos=" + str(pos))
                            poss[int(pos - 1)] += 1
                        if (compare):
                            dist = rf_distance.raxmlng_rf(
                                os.path.join(fam.get_gene_tree_dir(rep, family), tree_pick),
                                true_tree)
                            if (dist[0] == "abort"):
                                # print(seed, family, gt, '\n', dist)
                                writer.write(")\n")
                                continue
                            dist = dist[1]
                            # ", " if glob_pos=True
                            writer.write(", " * (glob_pos != None) + "dist=" + str(dist))
                            dists[int(dist * 100)][idx] += 1
                            if (dist <= 1):  # threshold
                                pick_avg_rel_dist[idx] = (pick_avg_rel_dist[idx] * pick_ctr[idx] + dist) / (pick_ctr[idx] + 1)
                                pick_ctr[idx] += 1
                        writer.write(")")
                    writer.write('\n')

    # metrics.update_dico(root_output, {"rank-distribution.rk" + str(i + 1) : poss[i] for i in range(len(poss))}, "misc")
    # metrics.update_dico(root_output, {"distance-distribution.dist" + str(i + 1): dists[i] for i in range(len(dists))}, "misc")
    with open(os.path.join(root_output, "metrics", tag + "_global__distributions.txt"), "w") as writer:
        writer.write("Rank distribution:")
        writer.write(str(poss))
        writer.write("\nDistance distribution (APro):")
        writer.write(str([dists[i][0] for i in range(len(dists))]))
        writer.write("\nDistance distribution (MAD):")
        writer.write(str([dists[i][1] for i in range(len(dists))]))
        writer.write("\nDistance distribution (All):")
        writer.write(str([dists[i][2] for i in range(len(dists))]))
    metrics.update_dico(root_output, {"pick_tree_avg_apro": pick_avg_rel_dist[0], "pick_tree_avg_mad": pick_avg_rel_dist[1], "pick_tree_avg_all": pick_avg_rel_dist[2]}, "misc")
    print("Avg pick distance (APro):", pick_avg_rel_dist[0])
    print("Avg pick distance (MAD):", pick_avg_rel_dist[1])
    print("Avg pick distance (All):", pick_avg_rel_dist[2])

def generax_likelihood_comp(root_output, replicates, tag, best_avg_tree, generax_run_dir):
    glob_picks_name = "_global__generax_picks.txt"
    stats_out = []
    for f in glob.glob(fam.get_metrics_dir(root_output) + "/" + tag + "*" + glob_picks_name):
        with open(f, 'r') as reader:
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
                try:
                    # stats = [raxmlLogLi, otherLogLi, sumLogLi, dup, loss]
                    stats_picked = launch_generax.eval(os.path.join(replicate, generax_run_dir, "results"), "family_" + family + ">" + tree_id)
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
    with open(os.path.join(root_output, "..", "logs", tag + "_likelihood_stats.txt"), 'w') as writer:
        for l in stats_out:
            writer.write(l)
