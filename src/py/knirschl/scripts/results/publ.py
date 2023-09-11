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

def distdistr2xy(vals):
    return ([i / 100 for i in range(len(vals))], vals)

def get_filename(tags, save):
    if not save:
        return None
    split = re.match(r"([A-Z]+)(.?[0-9.]+)?", tags[0])
    resultsdir = os.path.join("/home/fili/Desktop/2023/BA/code/output/benchmark_results/figures", ("bm", "single")[split[2] != None], tags[0])
    try:
        os.makedirs(resultsdir, exist_ok=True)
    except:
        pass
    return os.path.join(resultsdir, '_'.join(tags[1:]) + ".pdf")

# PLOTTER
def add2plot(plot, xy, ptype="scatter", labels=5 * [None], bottom=0, zoom=''):
    if (ptype == "scatter"):
        plt.make_scatter(plot, xy[0], xy[1], labels[0])
    elif (ptype == "bar"):
        plt.make_bar(plot, xy[0], xy[1], labels[0], labels[4], bottom)
    elif (ptype == "violin"):
        plt.make_violin(plot, xy[1], labels[0])
    else:
        plt.make_plot(xy[0], xy[1], label=labels[0])
    plt.set_titles(plot, title=labels[1], xAxis=labels[2], yAxis=labels[3])
    # TODO better zooming
    plt.cutoff(plot, xy[0], xy[1], zoom, THRESHOLD)
    
def plot_rrf(vals, tag, save):
    id = "rrf-dist"
    labels = ['', tag, "Skalierungswerte", "rRF-Distanz"]
    for zoom in ["auto", "hard"]:
        # plot all
        fig, ax = plt.subplots()
        for key in vals:
            labels[0] = key
            add2plot(ax, rrf2xy(vals, key), labels=labels, zoom=zoom)
        plt.display(fig, filename=get_filename([tag, "full", zoom, id], save))
        # plot single ba(+fm) variants
        for var in BA_VARIANTS:
            if (var in vals):
                #fig, ax = plt.subplots()
                labels[0:2] = [var, var + " (" + tag + ")"]
                add2plot(ax, rrf2xy(vals, var), labels=labels, zoom=zoom)
                plt.display(fig, makelegend=False, filename=get_filename([tag, var.lower(), zoom, id], save))

def plot_picks(vals, tag, save):
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
    plt.display(fig, filename=get_filename([tag, "generax-picks"], save))

def plot_distr(vals, tag, save):
    labels = [None, "Distanzen der von GeneRax ausgewählten Bäume (" + tag + ")", "rRF Distanz", "Anzahl an Bäumen", 0.01]
    # plot distances
    fig, ax = plt.subplots()
    for algo in ["APro","MAD", "All"]:
        labels[0] = algo
        ds = vals["distances (" + algo + ")"]
        add2plot(ax, distdistr2xy(ds[0]), ptype="bar", labels=labels)
        ax.axvline(ds[1] / 100, color='tab:green', label="avg (" + algo + ")")
    plt.display(fig, makegrid=False, filename=get_filename([tag, "dist-distr"], save))

# MAIN
def plot_single(dir, save=False):
    '''
    save == False -> show
    save == True -> don't show
    '''
    for f in os.listdir(dir):
        print(f)
        tag = f.split("_")[0]
        if "rf_distance" in f and "rel" in f:
            vals = reader.read_rrf(os.path.join(dir, f))
            plot_rrf(vals, tag, save)
        elif "generax_picks" in f:
            vals = reader.read_picks(os.path.join(dir, f))
            plot_picks(vals, tag, save)
        elif "distributions" in f:
            vals = reader.read_distr(os.path.join(dir, f))
            plot_distr(vals, tag, save)

def plot_bm(dir, save=False):
    '''
    save == False -> show
    save == True -> don't show
    '''
    # TODO plots (dist, runtimes, ...) of full benchmarks (e.g. [BASE, DUPLOS0.0, DUPLOS0.5, DUPLOS2.0, DUPLOS3.0])
    # TODO violin plots with error bars and not scatter
    # TODO collect them beforehand somehow? collect all in one dict and call plot_xx(..)

if (__name__ == "__main__"):
    plot_single("/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics", save=True)