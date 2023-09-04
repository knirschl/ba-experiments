import matplotlib.pyplot as plt
import numpy as np
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

COLORS = {
    GENERAX: "mediumblue",
    RAXML: "deepskyblue",
    FASTME: "darkviolet",
    build_ba_variant(BA, APRO): "lightgrey",
    build_ba_variant(BA, MAD): "darkgrey",
    build_ba_variant(BA, ALL): "dimgrey",
    build_ba_variant(BA_FASTME, APRO): "limegreen",
    build_ba_variant(BA_FASTME, MAD): "darkgreen",
    build_ba_variant(BA_FASTME, ALL): "yellowgreen",
    OTHER: "brown",
    None: None
}

MARKERS = {
    GENERAX: "o",
    RAXML: "x",
    FASTME: "d",
    build_ba_variant(BA, APRO): "o",
    build_ba_variant(BA, MAD): "x",
    build_ba_variant(BA, ALL): "d",
    build_ba_variant(BA_FASTME, APRO): "o",
    build_ba_variant(BA_FASTME, MAD): "x",
    build_ba_variant(BA_FASTME, ALL): "d",
    OTHER: "o",
    None: None
}

## HELPERS
def rrf2xy(vals, key):
    x = []
    y = []
    for stat in vals[key]:
        x.append(stat[2][0])
        y.append(stat[1])
    return (x, y)

def pickC2xy(vals, key, xs):
    x = []
    y = []
    for sc in xs:
        if (float(sc) <= 10): #zoom
            x.append(float(sc))
            if sc in vals[key]:
                y.append(int(vals[key][sc][0]))
            else:
                y.append(0.0)
    return (x, y)

def pickD2xy(vals, key):
    x = []
    y = []
    for sc in vals[key]:
        if (float(sc) <= 10): #zoom
            x.append(float(sc))
            y.append(float(vals[key][sc][1]))
    return (x, y)

def distdistr2xy(vals):
    return ([i / 100 for i in range(len(vals))], vals)

def get_filename(tags, save):
    if not save:
        return None
    return os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", '_'.join(tags) + ".pdf")

def add2plot(plot, xy, ptype="scatter", labels=5 * [None], bottom=0, zoom=''):
    if (ptype == "scatter"):
        plot.scatter(xy[0], xy[1], label=labels[0], c=COLORS[labels[0]], marker=MARKERS[labels[0]])
    elif (ptype == "bar"):
        plot.bar(xy[0], xy[1], label=labels[0], width=labels[4], color=COLORS[labels[0]], bottom=bottom)
    elif (ptype == "violin"):
        vp = plot.violinplot(xy[1], showmedians=True, showextrema=False, showmeans=False)
        vp['bodies'][0].set_edgecolor("k")
        vp['bodies'][0].set_facecolor(COLORS[labels[0]])
        vp['bodies'][0].set_alpha(0.75)
        for comp in ['cmedians']: #'cbars', 'cmins', 'cmaxes'
            vp[comp].set_edgecolor("lightgrey")
            vp[comp].set_alpha(1)
        plot.legend([vp['bodies'][0], vp['cmedians']], [labels[0], "median"])
    else:
        plot.plot(xy[0], xy[1], label=labels[0], color=COLORS[labels[0]], marker=MARKERS[labels[0]], markersize=2.5)
    plot.set(title=labels[1], xlabel=labels[2], ylabel=labels[3])
    # TODO better zooming
    if (zoom == "dyn"):
        plot.set_xlim([min(xy[0]) - 2, max(xy[0]) + 2])
        plot.set_ylim([min(xy[1]) - 0.01, max(xy[1]) + 0.01])
    elif (zoom == "hard"):
        plot.set_xlim([-0.2, 10.2])
        plot.set_ylim([0, 1])
    elif (zoom == "hardx"):
        plot.set_xlim([-0.2, 10.2])
    elif (zoom == "hardy"):
        plot.set_ylim([0, 1])

def display(figure, makegrid=True, makelegend=True, filename=None):
    plt.grid(makegrid)
    if makelegend:
        plt.legend()
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
    
# PLOTTER
def plot_rrf(rrf_file, tag, save):
    id = "rrf-dist"
    labels = ['', tag, "Skalierungswerte", "rRF-Distanz"]
    for zoom in ["auto", "hard"]:
        vals = read_rrf(rrf_file)
        # plot all
        fig, ax = plt.subplots()
        for key in vals:
            labels[0] = key
            add2plot(ax, rrf2xy(vals, key), labels=labels, zoom=zoom)
        display(fig, filename=get_filename([tag, "full", zoom, id], save))
        # plot single ba(+fm) variants
        for var in BA_VARIANTS:
            if (var in vals):
                #fig, ax = plt.subplots()
                labels[0:2] = [var, var + " (" + tag + ")"]
                add2plot(ax, rrf2xy(vals, var), labels=labels, zoom=zoom)
                display(fig, makelegend=False, filename=get_filename([tag, var.lower(), zoom, id], save))

def plot_picks(picks_file, tag, save):
    labels = ['', "GeneRax-Picks (" + tag + ")", "Skalierungswerte", "Häufigkeit", 0.05]
    vals = read_picks(picks_file)
    xs = set()
    for key in vals:
        for s in vals[key]:
            if (float(s) <= 10): #zoom
                xs.add(s)
    xy = {}
    for key in vals:
        xy[key] = list(zip(*sorted(zip(*pickC2xy(vals, key, xs)))))
    fig, ax = plt.subplots()
    bottoms = np.full(len(xs), 0)
    for key in xy:
        labels[0] = key
        # TODO zoom, correct labels
        add2plot(ax, xy[key], ptype="bar", labels=labels, bottom=bottoms)
        bottoms += np.array(xy[key][1], dtype='int64')
        """ labels[1] = None
        labels[3] = "rRF-Distanz"
        ax2 = ax.twinx()
        x, y = zip(*sorted(zip(*pickD2xy(vals, key))))
        add2plot(ax2, (x, y), ptype="plot", labels=labels) """
    display(fig, filename=get_filename([tag, "generax-picks"], save))

def plot_distr(distr_file, tag, save):
    labels = [None, "Distanzen der von GeneRax ausgewählten Bäume (" + tag + ")", "rRF Distanz", "Anzahl an Bäumen", 0.01]
    vals = read_distr(distr_file)
    # plot distances
    fig, ax = plt.subplots()
    ds = vals["distances"]
    add2plot(ax, distdistr2xy(ds[0]), ptype="bar", labels=labels)
    ax.axvline(ds[1] / 100, color='tab:green', label="avg")
    display(fig, makegrid=False, filename=get_filename([tag, "dist-distr"], save))

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
    # TODO plots with error bars and not scatter

if (__name__ == "__main__"):
    plot_matching_in_dir("/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics", save=True)