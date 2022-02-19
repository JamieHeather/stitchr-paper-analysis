# -*- coding: utf-8 -*-

"""
flowplotting.py

Plot the analysed flow cytometry data showing stitched-TCR+ Jurkat cell activavtion

"""


import os
import functions as fxn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import collections as coll
import math


__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.0'
__author__ = 'Jamie Heather'


def roundup(x, y):
    """
    :param x: a number to round
    :param y: a number to round up x to the nearest multiple of
    :return: x rounded up to nearest multiple of y
    """
    return int(math.ceil(x / y)) * y


# Set plotting parameters, read in data, filter out non-Jurkat controls
plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})
plot_dir = fxn.out_folder('flow-plotting')
best_dir = fxn.sub_out_folder(plot_dir, 'best-examples')
second_dir = fxn.sub_out_folder(plot_dir, 'second-best')
all_dat = pd.read_csv(fxn.data_dir + 'TCR-Jurkat-flow-data.tsv', sep='\t')
jurkats = ['C1-28', 'KFJ5', 'MAG-IC3', 'LC13', 'D2H']

j_dat = all_dat.loc[all_dat['Jurkat'].isin(jurkats)]
val = 'Jurkats | CD69+ CD62L- | %'  # The parameter that needs plotting

stimuli = ['0', '1', '10']  # omit 0.1 datapoints to save space, order left to right

pal = 'magma'
dat = coll.defaultdict()

# Pick the best plate to be shown per line
best_runs = {
    '2021-03-30',
    '2021-02-11',
    '2021-04-09',
    '2021-02-04',
    '2021-01-20'
}


# Scroll through data and plot
for tcr in jurkats:
    dat[tcr] = j_dat.loc[(j_dat['Jurkat'] == tcr) & (j_dat['HLA match'].isin(['Y', 'N']))]
    runs = list(set(dat[tcr]['Date']))
    for run in runs:
        # Plot the "best" runs into one directory, the second best into another
        if run in best_runs:
            specific_dir = best_dir
        else:
            specific_dir = second_dir

        out_str = tcr + '-' + run
        r_dat = dat[tcr].loc[dat[tcr]['Date'] == run]
        biggest = max(r_dat[val])
        if biggest < 40:
            roundto = 5
        else:
            roundto = 10
        upper = roundup(max(r_dat[val]), roundto)

        for suffix in ['.png']:
            fig = plt.figure(figsize=(2.5, 3.5))
            plt.rc('axes', linewidth=2)
            sns.boxplot(data=r_dat, x='HLA match', y=val, hue='Stimulus', palette=pal,
                        dodge=True, hue_order=stimuli, boxprops=dict(alpha=.2), order=['Y', 'N'])
            ax = sns.stripplot(data=r_dat, x='HLA match', y=val, hue='Stimulus', palette=pal, alpha=.7,
                          dodge=True, hue_order=stimuli, edgecolor='black', linewidth=1, size=10, order=['Y', 'N'])
            sns.despine(top=True, bottom=False, left=False, right=True)
            # plt.ylabel('% activated Jurkats', fontweight='bold')
            # plt.xlabel('HLA matched target', fontweight='bold')
            plt.ylabel('')
            plt.xlabel('')
            plt.ylim(0, upper)
            ax.get_legend().remove()
            # handles, labels = ax.get_legend_handles_labels()
            # l = plt.legend(handles[0:4], labels[0:4], bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
            plt.savefig(specific_dir + out_str + "-box-strip" + suffix, dpi=300, bbox_inches='tight')
            plt.close()

# And plot a mostly empty one, to crop out a legend
for suffix in ['.png']:
    fig = plt.figure(figsize=(2.5, 3.5))
    plt.rc('axes', linewidth=2)
    sns.boxplot(data=r_dat, x='HLA match', y=val, hue='Stimulus', palette=pal,
                dodge=True, hue_order=stimuli, boxprops=dict(alpha=.2), order=['Y', 'N'])
    sns.despine(top=True, bottom=True, left=True, right=True)
    plt.ylabel('')
    plt.xlabel('')
    plt.ylim(0, upper)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:3], labels[0:4], bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.savefig(plot_dir + "box-strip-legend" + suffix, dpi=300, bbox_inches='tight')
    plt.close()
