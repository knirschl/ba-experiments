import matplotlib.pyplot as plt
import re
import os

## CONSTANTS
def build_ba_variant(method, algorithm):
    return method + " (" + algorithm + ")"

# trees
GENERAX = "GeneRax"
RAXML = "RAxML"
FASTME = "FastME"
BA = "BA"
BA_FASTME = "BA+FM"
OTHER = "other"
# algos
APRO = "APro" 
MAD = "MAD"
ALL = "All"
ALGO_IDS = {
    'a': APRO,
    'm': MAD,
    '+': ALL
}
# subset of trees
BA_VARIANTS = [build_ba_variant(meth, a) for a in [APRO, MAD, ALL] for meth in [BA, BA_FASTME]]

""" COLORS = {
    GENERAX: "tab:blue",
    RAXML: "tab:purple",
    FASTME: "tab:orange",
    BA: "tab:red",
    BA_FASTME: "tab:green",
    OTHER: "tab:brown"
    } """

## HELPERS
def add_to_scatter(ax, vals, key, xlabel=None, ylabel=None, title=None):
    x = []
    y = []
    for stat in vals[key]:
        x.append(stat[2][0])
        y.append(stat[1])
    ax.scatter(x, y, label=key)#, c=COLORS[key])
    ax.set_title(title)
    ax.set(xlabel=xlabel, ylabel=ylabel)

def display(figure, makegrid=True, makelegend=True, filename=None):
    if makelegend:
        plt.legend()
    plt.grid(makegrid)
    if (filename):
        figure.savefig(filename, bbox_inches='tight')
    else:
        plt.show()
    plt.close(figure)

## PARSER
def read_rrf(rrf_file):
    '''
    Returns a dict containing every tree that has an entry in the file.
    Keys = {'generax', 'raxml', 'fastme', 'ba' + algo, 'ba+fastme' + algo (, 'other')}
    Values = [shortened tree name, rrf distance, [scale]
    Values for ba-trees = {}
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
    Returns ...
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
    
# PLOTTER
def plot_rrf(rrf_file, tag, save):
    vals = read_rrf(rrf_file)
    # plot all
    fig, ax = plt.subplots()
    for key in vals:
        add_to_scatter(ax, vals, key, "Skalierungswerte", "rRF-Distanz", tag)
    display(fig, filename=os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", tag + "_all_rf_distance-rel.pdf") if save else None)
    # plot single ba(+fm) variants
    for var in BA_VARIANTS:
        if (var in vals):
            fig, ax = plt.subplots()
            add_to_scatter(ax, vals, var, "Skalierungswerte", "rRF-Distanz", var + " (" + tag + ")")
            display(fig, makelegend=False, filename=os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", tag + "_" + var.lower() + "_rf_distance-rel.pdf") if save else None)

def plot_picks(picks_file, tag, save):
    vals = read_picks(picks_file)
    fig, ax = plt.subplots()
    for key in vals:
        x = []
        y = []
        for sc in vals[key]:
            if (float(sc) <= 10):
                x.append(float(sc))
                y.append(int(vals[key][sc][0]))
        # TODO stack
        ax.bar(x, y, label=key, width=0.05)#, c=COLORS[key])
        ax.set_title("GeneRax-Picks (" + tag + ")")
    display(fig, filename=os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", tag + "_generax_picks.pdf") if save else None)
        

def plot_distr(distr_file, tag, save):
    vals = read_distr(distr_file)
    # plot distances
    fig, ax = plt.subplots()
    ds = vals["distances"]
    ax.bar([i / 100 for i in range(len(ds[0]))], ds[0], width=0.01)
    ax.axvline(ds[1] / 100, color='tab:green', label="avg")
    ax.set_title("Distanzen der von GeneRax ausgewählten Bäume (" + tag + ")")
    ax.set(xlabel="rRF Distanz", ylabel="Anzahl an Bäumen")
    display(fig, makegrid=False, filename=os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", tag + "_distributions.pdf") if save else None)

# MAIN
def plot_matching_in_dir(dir, save=False):
    '''
    save == False -> show
    save == True -> don't show
    '''
    for f in os.listdir(dir):
        tag = f.split("_")[0]
        if "rf_distance" in f and "rel" in f:
            plot_rrf(os.path.join(dir, f), tag, save)
        elif "generax_picks" in f:
            plot_picks(os.path.join(dir, f), tag, save)
        elif "distributions" in f:
            plot_distr(os.path.join(dir, f), tag, save)
    # TODO plots (dist, runtimes, ...) of full benchmarks (e.g. [BASE, DUPLOS0.0, DUPLOS0.5, DUPLOS2.0, DUPLOS3.0])
    # TODO collect them beforehand somehow ?

if (__name__ == "__main__"):
    plot_matching_in_dir("/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics", save=False)