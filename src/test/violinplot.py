import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')


def inv(hex):
    return '#' + ''.join([hex(15 - int(c, 16))[2] for c in color[1:]])

def L(rgb):
    return (max(rgb) + min(rgb)) / 2

def luma(rgb):
    return 0.3 * rgb[0] + 0.59 * rgb[1] + 0.11 * rgb[2]

colors = mcolors.BASE_COLORS | mcolors.TABLEAU_COLORS | mcolors.CSS4_COLORS
Ls = {}
Lumas = {}
for c in colors:
    rgb = mcolors.to_rgb(c)
    Ls[c] = L(rgb)
    Lumas[c] = luma(rgb)
print("lightgrey:", Ls["lightgrey"], Lumas["lightgrey"])


# create test data
np.random.seed(19680801)
data = sorted(np.random.normal(0, 2, 100))
fig, ax = plt.subplots(nrows=13, ncols=13, figsize=(9, 4), sharey=True)

row = 0
column = 0
for color in colors:
    parts = ax[row][column].violinplot(data, showmedians=True, showextrema=False)
    line_color = ("lightgrey", "k")[abs(Lumas[color] - Lumas["lightgrey"]) < 0.4]
    parts['bodies'][0].set_facecolor(color)
    parts['bodies'][0].set_edgecolor('k')
    parts['bodies'][0].set_alpha(1)

    #parts['cmins'].set_edgecolor(line_color)
    #parts['cmaxes'].set_edgecolor(line_color)
    #parts['cbars'].set_edgecolor(line_color)
    parts['cmedians'].set_edgecolor(line_color)

    if row < 12:
        row += 1
    else:
        row = 0
        column += 1 
#plt.subplots_adjust(bottom=0.15, wspace=0.05)
plt.show()