import numpy as np
import re
import os
import reader
from constants import *
import plots.plotter as plt

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
        if (float(sc) <= THRESHOLD): #zoom
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
        if (float(sc) <= THRESHOLD): #zoom
            x.append(float(sc))
            y.append(float(vals[key][sc][1]))
    return (x, y)

def extract_xy(vals, tool, idx):
    xs = []
    ys = []
    for x in vals:
        xs.append(float(x))
        ys.append(sum(vals[x][tool][idx]) / ([1, len(vals[x][tool][idx])][idx == 1]))
    return (xs, ys)

def distdistr2xy(vals):
    return ([i / 100 for i in range(len(vals))], vals)

def get_filename(tags, save):
    if not save:
        return None
    split = re.match(r"([A-Z]+)(.?[0-9.]+)?", tags[1])
    resultsdir = os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", tags[0], ("bm", "single")[split[2] != None], tags[1])
    try:
        os.makedirs(resultsdir, exist_ok=True)
    except:
        pass
    return os.path.join(resultsdir, '_'.join(tags[2:]) + ".pdf")

# PLOTTER
def add2plot(plot, xy, ptype="scatter", labels=5 * [None], bottom=0, zoom=''):
    if (ptype == "scatter"):
        plt.make_scatter(plot, xy[0], xy[1], labels[0])
    elif (ptype == "bar"):
        plt.make_bar(plot, xy[0], xy[1], labels[0], labels[4], bottom)
    elif (ptype == "violin"):
        plt.make_violin(plot, xy[1], labels[0])
    else:
        xy = tuple((zip(*sorted(zip(*xy)))))
        plt.make_plot(plot, xy[0], xy[1], label=labels[0], logscale=(labels[1] == "BRALEN" and "Skalier" not in labels[2]))
    plt.set_titles(plot, title=labels[1], xAxis=labels[2], yAxis=labels[3])
    plt.cutoff(plot, xy[0], xy[1], zoom, THRESHOLD)
    
def plot_rrf_single(vals, setup, tag, save):
    id = "rrf-dist"
    labels = ['', tag, "Skalierungswerte", "rRF-Distanz"]
    for zoom in ["auto", "hard"]:
        # plot all
        fig, ax = plt.subplots()
        for key in vals:
            labels[0] = key
            add2plot(ax, rrf2xy(vals, key), labels=labels, zoom=zoom)
        plt.display(fig, filename=get_filename([setup, tag, "full", zoom, id], save))
        # plot single ba(+fm) variants
        for var in BA_VARIANTS:
            if (var in vals):
                #fig, ax = plt.subplots()
                labels[0:2] = [var, var + " (" + tag + ")"]
                add2plot(ax, rrf2xy(vals, var), labels=labels, zoom=zoom)
                plt.display(fig, makelegend=False, filename=get_filename([setup, tag, var.lower(), zoom, id], save))

def plot_rrf(bm, setup, tag, save):
    id = "rrf-distance"
    labels = ['', tag, "Varianten", "rRF-Distanz"]
    fig, ax = plt.subplots()
    ticks = set()
    for tool in bm:
        labels[0] = tool
        x = []
        y = []
        for var in bm[tool]:
            ticks.add(float(var))
            x.append(float(var))
            y.append(bm[tool][var])
        add2plot(ax, (x, y), labels=labels, ptype="plot")
    ax.set_xticks(sorted(list(ticks)))
    plt.display(fig, filename=get_filename([setup, tag, tag.lower(), id], save))

def plot_rt(bm, setup, tag, save):
    id = "runtimes"
    labels = [None, tag, "Varianten", "Laufzeit in Sekunden"]
    fig, ax = plt.subplots()
    for tool in bm:
        labels[0] = tool
        x = []
        y = []
        for var in bm[tool]:
            x.append(float(var))
            y.append(bm[tool][var])
        add2plot(ax, (x, y), labels=labels, ptype="plot")
    plt.display(fig, filename=get_filename([setup, tag, tag.lower(), id], save))

def plot_pick(bm, setup, tag, save):
    vals = {}
    xs = set()
    # pre-compute
    for tool in bm:
        for var in bm[tool]:
            btv = bm[tool][var]
            for scale in btv:
                if (float(scale) > THRESHOLD): # "zoom"
                    continue
                if scale not in vals:
                    xs.add(scale)
                    vals[scale] = {}
                if tool not in vals[scale]:
                    vals[scale][tool] = [[], []] # [occurences, avg distance]
                vals[scale][tool][0].append(btv[scale][0])
                vals[scale][tool][1].append(btv[scale][1])
    # plot occurences
    id = "picks_occurences"
    labels = [None, tag, "Skalierungswerte", "Häufigkeit", 0.05]
    fig, ax = plt.subplots()
    bottoms = np.full(len(xs), 0)
    for tool in bm:
        labels[0] = tool
        xy = extract_xy(vals, tool, 0)
        add2plot(ax, xy, labels=labels, ptype="bar", bottom=bottoms)
        bottoms += np.array(xy[1], dtype='int64')
    plt.display(fig, filename=get_filename([setup, tag, tag.lower(), id], save))
    # plot average distances
    id = "picks_avg_dist"
    labels = [None, tag, "Skalierungswerte", "rRF-Distanz"]
    fig, ax = plt.subplots()
    for tool in bm:
        labels[0] = tool
        add2plot(ax, extract_xy(vals, tool, 1), labels=labels, ptype="plot")
    plt.display(fig, filename=get_filename([setup, tag, tag.lower(), id], save))

def plot_picks_single(vals, setup, tag, save):
    # TODO still working after modifying generax picking ??
    
    labels = ['', "GeneRax-Picks (" + tag + ")", "Skalierungswerte", "Häufigkeit", 0.05]
    xs = set()
    for key in vals:
        for s in vals[key]:
            if (float(s) <= THRESHOLD): #zoom
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
    plt.display(fig, filename=get_filename([setup, tag, "generax-picks"], save))

def plot_distr_single(vals, setup, tag, save):
    labels = [None, "Distanzen der von GeneRax ausgewählten Bäume (" + tag + ")", "rRF Distanz", "Anzahl an Bäumen", 0.01]
    # plot distances
    fig, ax = plt.subplots()
    for algo in ["APro","MAD", "All"]:
        labels[0] = algo
        ds = vals["distances (" + algo.lower() + ")"]
        add2plot(ax, distdistr2xy(ds[0]), ptype="bar", labels=labels)
        ax.axvline(ds[1] / 100, color='tab:green', label="avg (" + algo + ")")
    plt.display(fig, makegrid=False, filename=get_filename([setup, tag, "dist-distr"], save))

# MAIN
def plot_single(dir, save=False):
    '''
    save == False -> show
    save == True -> don't show
    '''
    for test_setup in os.listdir(dir):
        for f in os.listdir(os.path.join(dir, test_setup)):
            print(f)
            tag = f.split("_")[0]
            if "rf_distance" in f and "rel" in f:
                vals = reader.read_rrf(os.path.join(dir, test_setup, f))
                plot_rrf_single(vals, test_setup, tag, save)
            elif "generax_picks" in f:
                vals = reader.read_picks(os.path.join(dir, test_setup, f))
                plot_picks_single(vals, test_setup, tag, save)
            elif "distributions" in f:
                vals = reader.read_distr(os.path.join(dir, test_setup, f))
                plot_distr_single(vals, test_setup, tag, save)

def plot_bm(dir, save=False):
    '''
    save == False -> show
    save == True -> don't show
    '''
    for test_setup in os.listdir(dir):
        # benchmark : {variation : vals}
        ds_rrf = {k: {} for k in DATASETS} # rf distance
        ds_rt = {k: {} for k in DATASETS} # runtimes
        ds_pick = {k: {} for k in DATASETS} # picks
        for f in os.listdir(os.path.join(dir, test_setup)):
            rem = re.match(r"([A-Z]+)(-?[0-9.]+)_global__([a-z_-]+).txt", f)
            ds_name = rem[1]
            ds_var = rem[2]
            metric = rem[3]
            if metric == "rf_distance_avg-rel":
                ds_results = ds_rrf
                read = reader.read_arf_pick(os.path.join(dir, test_setup, f), os.path.join(dir, test_setup, f.replace("rf_distance_avg-rel", "generax_picks")))
            elif metric == "runtimes_avg":
                ds_results = ds_rt
                read = reader.read_rt(os.path.join(dir, test_setup, f), scaling=test_setup[:test_setup.index('x')]) # set scaling to number of tested methods (2x?? -> 2)
            elif metric == "generax_picks":
                ds_results = ds_pick
                read = reader.read_picks(os.path.join(dir, test_setup, f))
            else:
                continue
            if (ds_name != "BASE"):
                ds_results[ds_name][ds_var] = read
            else:
                for k in DATASETS:
                    ds_results[k][DATASETS[k]] = read
        # plot
        for benchmark in ds_rrf:
            plot_rrf(reader.map_bm(ds_rrf[benchmark]), test_setup, benchmark, save)
        for benchmark in ds_rt:
            plot_rt(reader.map_bm(ds_rt[benchmark]), test_setup, benchmark, save)
        for benchmark in ds_pick:
            #print(benchmark, reader.map_bm(ds_pick[benchmark]))
            plot_pick(reader.map_bm(ds_pick[benchmark]), test_setup, benchmark, save)

if (__name__ == "__main__"):
    #plot_single("/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics", save=True) # plots look broken?
    plot_bm("/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics", save=True)