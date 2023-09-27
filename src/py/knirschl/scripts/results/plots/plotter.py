from plots.settings import *
import matplotlib.pyplot as plt

def subplots(nrows=1, ncols=1, *, sharex=False, sharey=False, squeeze=True, width_ratios=None, height_ratios=None, subplot_kw=None, gridspec_kw=None, **fig_kw):
    return plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey, squeeze=squeeze, width_ratios=width_ratios, height_ratios=height_ratios, subplot_kw=subplot_kw, gridspec_kw=gridspec_kw, **fig_kw)

def make_plot(plot, x, y, label, logscale=False):
    raw_label = label.split('$')[0]
    plot.plot(x, y, label=label, color=COLORS[raw_label], marker=MARKERS[raw_label])#, markersize=2.5)
    plt.xscale(("linear", "log")[logscale])

def make_scatter(plot, x, y, label):
    raw_label = label.split('$')[0]
    plot.scatter(x, y, label=label, c=COLORS[raw_label], marker=MARKERS[raw_label])

def make_bar(plot, x, y, label, width, bottom):
    raw_label = label.split('$')[0]
    plot.bar(x, y, label=label, color=COLORS[raw_label], width=width, bottom=bottom)

def make_violin(plot, y, label):
    vp = plot.violinplot(y, showmedians=True, showextrema=False, showmeans=False)
    vp['bodies'][0].set_edgecolor("k")
    vp['bodies'][0].set_facecolor(COLORS[label])
    vp['bodies'][0].set_alpha(0.75)
    for comp in ['cmedians']: #'cbars', 'cmins', 'cmaxes'
        vp[comp].set_edgecolor("lightgrey")
        vp[comp].set_alpha(1)
    plot.legend([vp['bodies'][0], vp['cmedians']], [label, "median"])

def display(figure, makegrid=True, makelegend=True, filename=None):
    plt.grid(makegrid)
    plt.ylim(ymin=0)
    if makelegend:
        plt.legend(framealpha=1)#, loc='lower right')
    if (filename):
        figure.savefig(filename, bbox_inches='tight')
    else:
        plt.show()
    plt.close(figure)

def set_titles(plot, title=None, xAxis=None, yAxis=None):
    plot.set(title=title, xlabel=xAxis, ylabel=yAxis)

def cutoff(plot, x, y, zoom, thresh):
    if (zoom == "dyn"):
        plot.set_xlim([min(x) - 2, max(x) + 2])
        plot.set_ylim([min(y) - 0.01, max(y) + 0.01])
    elif (zoom == "hard"):
        plot.set_xlim([-0.2, thresh + 0.2])
        plot.set_ylim([0, 1])
    elif (zoom == "hardx"):
        plot.set_xlim([-0.2, thresh + 0.2])
    elif (zoom == "hardy"):
        plot.set_ylim([0, 1])