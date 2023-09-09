import re
from constants import *

def read_rrf(rrf_file):
    '''
    Returns a dict containing every tree that has an entry in the file.
    Keys = {'generax', 'raxml', 'fastme', 'ba' + algo, 'ba+fastme' + algo (, 'other')}
    Values = [shortened tree name, rrf distance, [scale (, fastme matrix model, starting tree)]
    '''
    trees = {
        GENERAX: [],
        RAXML: [],
        FASTME: []
    }
    for var in BA_VARIANTS:
        trees[var] = []
    with open(rrf_file, 'r') as reader:
        for line in reader.readlines():
            current_tree, current_dist = line.split(" : ")
            current_tree = current_tree.replace(".genetree", '').replace(".newick", '') # fastme trees doesn't contain '.geneTree'
            current_dist = float(current_dist)
            if (GENERAX.lower() in current_tree):
                trees[GENERAX].append((current_tree, current_dist, [0.0]))
            elif (RAXML.lower() in current_tree):
                trees[RAXML].append((current_tree, current_dist, [0.0]))
            elif (BA.lower() in current_tree):
                groups = re.match(r"ba\.(\S+)_(\S+)([am+])\.(?:fastme\.)?([0-9.]+)\S+", current_tree)
                fastme_model = groups[1]
                starting_tree = groups[2]
                reroot_algo = groups[3]
                scale = float(groups[4])
                if (FASTME.lower() in current_tree):
                    var = BA_FASTME
                else:
                    var = BA
                trees[build_ba_variant(var, ALGO_IDS[reroot_algo])].append((current_tree, current_dist, [scale, fastme_model, starting_tree]))                
            elif (FASTME.lower() in current_tree):
                trees[FASTME].append((current_tree, current_dist, [0.0]))
            else:
                if OTHER not in trees:
                    trees[OTHER] = []
                trees[OTHER].append((current_tree, current_dist, [0.0]))
    trees = {key: val for key, val in trees.items() if val != []}
    #print(trees)
    return trees

def read_picks(picks_file):
    '''
    Returns a dict containing dicts per used scale.
    Keys = {'ba' + algo, 'ba+fastme' + algo}
    Values = {scaling values}
    Sec. Valuies = [no. occurences, average distance, [dists], [trees]]
    '''
    trees = {}
    for var in BA_VARIANTS:
        trees[var] = {}
    with open(picks_file, 'r') as reader:
        for line in reader.readlines():
            if re.match(r"[0-9]+", line):
                # seed
                continue
            _, current_tree, stats = line.split("  ")
            current_dist = float(re.match(r".*dist=([0-9.]+).*", stats)[1])
            groups = re.match(r"ba\.(\S+)_(\S+)([am+])\.(?:fastme\.)?([0-9.]+)\S+", current_tree)
            reroot_algo = groups[3]
            scale = groups[4]
            if (FASTME.lower() in current_tree):
                var = BA_FASTME
            else:
                var = BA
            if scale not in trees[build_ba_variant(var, ALGO_IDS[reroot_algo])]:
                trees[build_ba_variant(var, ALGO_IDS[reroot_algo])][scale] = [0, -1, [], []]
            trees[build_ba_variant(var, ALGO_IDS[reroot_algo])][scale][0] += 1
            trees[build_ba_variant(var, ALGO_IDS[reroot_algo])][scale][2].append(current_dist)
            trees[build_ba_variant(var, ALGO_IDS[reroot_algo])][scale][3].append(current_tree)
    trees = {key: val for key, val in trees.items() if val != {}}
    for key in trees:
        for scale in trees[key]:
            curr = trees[key][scale]
            curr[1] = sum(curr[2]) / len(curr[2])
    #print(trees)
    return trees

def read_distr(distr_file):
    '''
    Returns a dict containing every distribution in the file.
    Values = (distribution, average (weighted by index))
    '''
    distributions = {}
    with open(distr_file, 'r') as reader:
        for line in reader.readlines():
            name, ls = line.split(" distribution:")
            name = name.lower() + 's'
            ds = []
            for v in ls.replace('[', '').replace(']', '').split(", "):
                ds.append(int(v))
            distributions[name] = (ds, sum([ds[i] * i for i in range(len(ds))]) / sum(ds))
    #print(distributions)
    return distributions