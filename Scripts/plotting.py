# -*- coding: utf-8 -*-

"""
plotting.py

Run the actual plotting that underlies the Stitchr manuscript.

"""


import os
import collections as coll
import functions as fxn
import Levenshtein as lev
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")


__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.2'
__author__ = 'Jamie Heather'


def fix_line(line_in):
    """
    :param line_in: Str line of TCR file to be sanitised
    :return: line_in as list split on tabs, but with quotation and return characters removed
    """
    return line_in.replace('\"','').rstrip().split('\t')


def fix_post_seq(post_seq, pre_seq_len):
    """
    :param post_seq: TCR nt produced by stitchr
    :param pre_seq_len: length of inter-tag sequence from autoDCR
    :return: Sequence with constant region and artificial 5' ATG removed, running up to length of autoDCR seq
    """
    for c in [trac, trbc1, trbc2]:
        post_seq = post_seq.replace(c, "")
    return post_seq[3:pre_seq_len+3]


def get_stitchr_code(input_tcr):
    """
    Function for manual QC checks of TCRs from different data formats generated during analysis
    :param input_tcr: a TCR from one of the data structures they're stored in throughout this script
    :return: Nothing, but it prints text which can be supplied to run stitchr in the terminal
    """
    if isinstance(input_tcr, list):
        print(' '.join(['python3 stitchr.py -v', input_tcr[1],
                        '-j', input_tcr[2],
                        '-cdr3', input_tcr[3]]))
    elif isinstance(input_tcr, pd.Series):
        print(' '.join(['python3 stitchr.py -v', input_tcr['v_call'],
                        '-j', input_tcr['j_call'],
                        '-cdr3', input_tcr['junction']]))
    else:
        raise IOError("Type not covered.")


def bin_errors(num_errors):
    """
    :param num_errors: Integer of number of errors in a given read
    :return: one of several categories to group different approximate error ranges together
    """
    if not isinstance(num_errors, int):
        raise ValueError("Not an integer.")
    elif num_errors < 0:
        raise ValueError("Cannot process negative integers.")
    elif num_errors == 0:
        return '0'
    elif num_errors <= 2:
        return '1-2'
    elif num_errors <= 5:
        return '3-5'
    else:
        return '6+'


def plot_multi_modal(dir_to_thimble, dir_to_proc, file_prefix, out_dir_nam, comparison_file, whole_seq_field):
    """
    :param dir_to_thimble: Path to directory containing the stitched files output by Thimble
    :param dir_to_proc: Path to directory containing the pre-Stitchr/Thimble processed TCR data
    :param file_prefix: Str of input file prefix
    :param out_dir_nam: Str name to distinguish output sub-folder
    :param comparison_file: Str name of file in dir_to_proc containing the pre-stitched TCR information
    :param whole_seq_field: Str olumn ID in the comparison_file that contains the actual raw input TCR sequence
    :return: all_dat pandas dataframe, containing the TCR and error info of all TCRs processed in this dataset
    """

    # The major plotting task: comparing the error profiles of the different stitching options against full length seqs
    out_dir = fxn.sub_out_folder(base_out_dir, out_dir_nam)
    stitched_files = [x for x in os.listdir(dir_to_thimble) if 'timing_data' not in x
                      and '.tsv' in x and file_prefix in x and '~' not in x]

    stitched_files.sort()

    # dictionaries of dataframes for each analysis
    errs = coll.defaultdict()
    all_dat = []
    all_rel_pos = []
    all_rel_pos_aa = []
    indels = []

    pre_stitched_file = dir_to_proc + [x for x in os.listdir(dir_to_proc) if comparison_file in x][0]

    print("Calculating " + file_prefix + " dataset error profiles relative to 'true' reads")

    for stitch_file in stitched_files:

        nam = stitch_file.split('.')[0].\
            replace(file_prefix + '_NT-SL-10', 'SL010').\
            replace(file_prefix + '_NT-SL-200', 'SL200').\
            replace(file_prefix + '_NT-SL-20', 'SL020').\
            replace(file_prefix + '_NT-SL-30', 'SL030').\
            replace(file_prefix + '_AA', 'AA').\
            replace(file_prefix + '_NT', 'NT')

        print('Name:\t' + nam)

        errors_df = []
        hdists = coll.Counter()
        ldists = coll.Counter()
        ratios = []

        # Loop through both files simultaneously, reading in appropriate data
        with fxn.opener(pre_stitched_file, 'r') as pre_in, fxn.opener(dir_to_thimble + stitch_file, 'r') as post_in:
            line_count = 0
            for pre_line, post_line in zip(pre_in, post_in):

                # Get/check basic line information
                pre = fix_line(pre_line)
                post = fix_line(post_line)

                if line_count == 0:
                    pre_headers = pre
                    post_headers = post

                else:
                    # Check we're in sync
                    if pre[0] != post[0]:
                        print(pre, post)
                        print(pre_stitched_file, dir_to_thimble + stitch_file)
                        raise IOError("Files not in sync! See line number " + str(line_count + 1))

                    # Then pull out the relevant TCR info
                    tcr_id = pre[0]
                    if 'sequence' in pre_headers:
                        tcr_seq = pre[pre_headers.index('sequence')]
                    else:
                        tcr_seq = pre[pre_headers.index('inter_tag_seq')]
                    v = pre[pre_headers.index('v_call')]
                    j = pre[pre_headers.index('j_call')]
                    jxn = pre[pre_headers.index('junction')]

                    donor_line = pre[7].replace('\'', '').replace('[', '').\
                        replace(']', '').replace(' ', '').split(',')
                    donors = [x.split('-')[0] for x in donor_line]
                    donor_freqs = [x.split(':')[-1] for x in donor_line]

                    tcr_chain = pre[pre_headers.index('j_call')][:3]
                    whole_raw_seq = pre[pre_headers.index(whole_seq_field)].upper()
                    stitched_seq = post[post_headers.index(tcr_chain + '_nt')]
                    errors = post[post_headers.index('Warnings/Errors')]

                    if "Cannot stitch a sequence" in errors:
                        errors_df.append(pre + [errors])
                        continue
                    elif '*' in fxn.translate_nt(whole_raw_seq):
                        errors_df.append(pre + [errors + ' AUTODCR FAILURE!'])
                        continue

                    # Fix the stitched TCR, to get equivalent sequences as Decombined
                    matched_stitched = fix_post_seq(stitched_seq, len(whole_raw_seq))
                    if len(whole_raw_seq) == len(matched_stitched):
                        hdist = lev.hamming(whole_raw_seq, matched_stitched)
                    else:
                        hdist = np.nan

                    ldist = lev.distance(whole_raw_seq, matched_stitched)
                    ratio = lev.ratio(whole_raw_seq, matched_stitched)
                    hdists[hdist] += 1
                    ldists[ldist] += 1
                    ratios.append(ratio)

                    # And then do the same for amino acids
                    matched_stitched_aa = fxn.translate_nt(matched_stitched)
                    whole_raw_seq_aa = fxn.translate_nt(whole_raw_seq)
                    ldist_aa = lev.distance(whole_raw_seq_aa, matched_stitched_aa)
                    hdist_aa = lev.distance(whole_raw_seq_aa, matched_stitched_aa)

                    # We want to gather the number of nucleotides/AA, for residue-level metrics
                    num_nt = len(whole_raw_seq)
                    num_nt_right = num_nt - hdist
                    num_aa = len(whole_raw_seq_aa)
                    num_aa_right = num_aa - hdist_aa

                    # If there's an error, calculate its relative position
                    if ldist > 0:
                        if len(whole_raw_seq) == len(matched_stitched):
                            # Both for nucleotides...
                            for nt in range(len(whole_raw_seq)):
                                if whole_raw_seq[nt] != matched_stitched[nt]:
                                    rel_pos = nt / len(whole_raw_seq)
                                    all_rel_pos.append([nam, tcr_id, tcr_chain, rel_pos, ldist])

                            # ... and for amino acids
                            for aa in range(len(whole_raw_seq_aa)):
                                if whole_raw_seq_aa[aa] != matched_stitched_aa[aa]:
                                    rel_pos = aa / len(whole_raw_seq_aa)
                                    all_rel_pos_aa.append([nam, tcr_id, tcr_chain, rel_pos, ldist_aa])

                        else:
                            indels.append([nam, tcr_id, tcr_chain, 1])

                    # and save details of pan-mode data for cross-comparisons
                    all_dat.append([nam, tcr_id, tcr_chain, v, j, jxn, tcr_seq, donors, donor_freqs,
                                    ldist, hdist, ratio, num_nt, num_nt_right,
                                    bin_errors(hdist), bin_errors(hdist_aa),
                                    ldist_aa, hdist_aa, num_aa, num_aa_right])

                line_count += 1

        if errors_df:
            errors = fxn.list_to_df(errors_df, pre_headers + ['errors'], False)
        else:
            errors = fxn.list_to_df([' ' for x in range(len(pre_headers) + 1)], pre_headers + ['errors'], False)

        errs[nam] = errors

    out_headers = ['Run', 'TCR_ID', 'Chain', 'V', 'J', 'Junction', 'Sequence', 'Donors', 'Donor_Freqs',
                   'Edit_Distance', 'Hamming_Distance', 'Similarity_Ratio',
                   'Number_NT', 'Number_NT_Correct', 'Number_Errors_NT', 'Number_Errors_AA',
                   'Edit_Distance_AA', 'Hamming_Distance_AA', 'Number_AA', 'Number_AA_Correct']

    all_dat = fxn.list_to_df(all_dat, out_headers, False)

    all_dat.to_csv(out_dir + 'all_dat.tsv.gz', sep='\t', compression='gzip')

    if all_rel_pos:
        all_rel_pos = fxn.list_to_df(all_rel_pos,
                                     ['Run', 'TCR_ID', 'Chain', 'Relative_Position', 'Edit_Distance'], False)
        all_rel_pos = seamless_rename(all_rel_pos)

    if all_rel_pos_aa:
        all_rel_pos_aa = fxn.list_to_df(all_rel_pos_aa,
                                        ['Run', 'TCR_ID', 'Chain', 'Relative_Position', 'Edit_Distance'], False)
        all_rel_pos_aa = seamless_rename(all_rel_pos_aa)

    if indels:
        indels = fxn.list_to_df(indels, ['Run', 'TCR_ID', 'Chain', 'Count'], False)  # Note all counts are '1'
    else:
        print("No insertion/deletion errors detected!")

    # Plot nucleotide-level metrics
    runs = list(set(all_dat['Run']))
    runs.sort()

    res_dat = []
    for r in runs:
        tmp_dat = all_dat.loc[all_dat['Run'] == r]
        for typ in ['NT', 'AA']:
            percent_res_right = (sum(tmp_dat['Number_' + typ + '_Correct']) / sum(tmp_dat['Number_' + typ])) * 100
            edit_str = 'Edit_Distance'
            if typ == 'AA':
                edit_str += '_AA'
            percent_tcr_perfect = (len(tmp_dat.loc[tmp_dat[edit_str] == 0]) / len(tmp_dat)) * 100
            res_dat.append([r, typ, percent_res_right, percent_tcr_perfect])

    res_dat = fxn.list_to_df(res_dat, ['Run', 'Type', 'Percent Residues Identical', 'Percent TCRs Perfect'], False)

    # Put all the carriage returns back in the run names, for plotting
    all_dat = seamless_rename(all_dat)
    res_dat = seamless_rename(res_dat)

    for ext in fxn.exts:

        for sz in fxn.sizes:
            sns.displot(data=all_dat, kind='kde', x='Number_NT', hue='Chain', hue_order=chains, height=sz)
            plt.xlim(250, 400)
            plt.ylabel('Density', fontweight='bold')
            plt.xlabel('Length of variable domain', fontweight='bold')
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'variable-domain-length-kde.' +
                        str(sz).replace('.', '-') + '.' + ext, dpi=300)  #, bbox_inches='tight')
            plt.close('all')

            fig = plt.figure(figsize=(sz, sz))
            sns.countplot(data=all_dat, x='Run', hue='Number_Errors_NT', hue_order=error_bins,
                          order=cdr3_order, palette=sns.color_palette('bright')[2:])
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            sns.despine(top=True, bottom=False, left=False, right=True)
            plt.xlabel("")
            plt.ylabel('Number of TCRs', fontweight='bold')
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'number-errors-NT-countplot.' +
                        str(sz).replace('.', '-') + '.' + ext, dpi=300, bbox_inches='tight')
            plt.close('all')

            fig = plt.figure(figsize=(sz, sz))
            sns.countplot(data=all_dat, x='Run', hue='Number_Errors_AA', hue_order=error_bins,
                          order=cdr3_order, palette=sns.color_palette('bright')[2:])
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            sns.despine(top=True, bottom=False, left=False, right=True)
            plt.xlabel("")
            plt.ylabel('Number of TCRs', fontweight='bold')
            plt.savefig(out_dir + 'number-errors-AA-countplot.' +
                        str(sz).replace('.', '-') + '.' + ext, dpi=300, bbox_inches='tight')
            plt.close('all')

            sns.catplot(data=res_dat, kind='bar', x='Run', hue='Type', y='Percent Residues Identical', height=sz,
                        order=cdr3_order_plotting_full, palette=cdr3_palette)
            plt.ylim([98.5, 100])
            plt.axhline(y=100, alpha=.5, color='black', ls='--', lw=3)
            plt.xlabel("")
            plt.ylabel('Identical residues (%)', fontweight='bold')
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'percent-residues-identical.' +
                        str(sz).replace('.', '-') + '.' + ext, dpi=300)  #, bbox_inches='tight')
            plt.close('all')

            fig = plt.figure(figsize=(sz, sz))
            g = sns.catplot(data=res_dat, kind='bar', x='Run', hue='Type', y='Percent TCRs Perfect', height=sz,
                            hue_order=res_order, order=cdr3_order_plotting_full, palette=cdr3_palette)
            plt.ylim([0, 100])
            # plt.axhline(y=100, alpha=.5, color='black', ls='--', lw=3)
            plt.xlabel("")
            plt.ylabel('Perfectly matching TCRs (%)', fontweight='bold')
            plt.rc('axes', linewidth=2)

            offsets = [-0.18, 0.23]
            ax = g.axes[0, 0]

            # Plot text values of percentages
            for x in range(len(cdr3_order_plotting_full)):
                row = res_dat.loc[res_dat['Run'] == cdr3_order_plotting_full[x]]
                for res_pos in range(len(res_order)):
                    res = res_order[res_pos]
                    y_val = row.loc[row['Type'] == res]['Percent TCRs Perfect'][0]
                    if y_val < 20:
                        plot_y = y_val + 5
                        txt_col = cdr3_palette[res_pos]
                    else:
                        plot_y = y_val - 15
                        txt_col = 'white'
                    if y_val == 100:
                        txt_y = str(100)
                    else:
                        txt_y = round(y_val, 1)
                    ax.text(x + offsets[res_pos], plot_y, txt_y, color=txt_col,
                            ha="center", fontsize=14, rotation='vertical', horizontalalignment='center')

            plt.savefig(out_dir + 'percent-TCRs-perfect-' + str(sz).replace('.', '-') + '.' + ext,
                        dpi=300)  #, bbox_inches='tight')
            plt.close('all')

            nbins = 111         # 111 = about 1/3 of the average TCR length, so 1 bar = ~ 1 codon (easier to see)
            bw = 1/nbins
            if len(cdr3_order) > 1:
                used_runs =  [x for x in cdr3_order_plotting if x in list(set(all_rel_pos['Run']))]

                fig = plt.figure(figsize=(sz*4, sz))
                g = sns.displot(data=all_rel_pos, x='Relative_Position', kind='hist', col='Run', hue='Chain', height=sz,
                                rug=True, rug_kws={"alpha": 0.3}, binwidth=bw, col_order=used_runs, hue_order=chains)
                g.set_xlabels('Relative Position in TCR', fontweight='bold')
                g.set_ylabels('Number Errors', fontweight='bold')
                axes = g.axes.flatten()
                for i in range(len(used_runs)):
                    axes[i].set_title(used_runs[i], fontweight='bold', verticalalignment='bottom')
                    axes[i].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
                    axes[i].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

                ax = g.axes[0, 0]
                max_val = np.max(ax.get_ylim()[1])
                plt.gca().set_ylim(bottom=-(max_val / 30))
                plt.xlim(0, 1)
                plt.rc('axes', linewidth=2)
                sns.despine(top=False, bottom=False, left=False, right=False)
                plt.savefig(out_dir + 'relative-NT-error-hist-' + str(sz).replace('.', '-') + '.' + ext,
                            dpi=300)
                plt.close('all')

                bottom_hundredth = max_val / 100

                fig = plt.figure(figsize=(sz*4, sz))
                g = sns.displot(data=all_rel_pos, x='Relative_Position', kind='hist', col='Run', hue='Chain', height=sz,
                                rug=True, rug_kws={"alpha": 0.3}, binwidth=bw, col_order=used_runs, hue_order=chains)
                g.set_xlabels('Relative Position in TCR', fontweight='bold')
                g.set_ylabels('Number Errors', fontweight='bold')
                axes = g.axes.flatten()
                for i in range(len(used_runs)):
                    axes[i].set_title(used_runs[i], fontweight='bold', verticalalignment='bottom')
                    axes[i].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
                    axes[i].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

                ax = g.axes[0, 0]
                plt.gca().set_ylim(bottom=-(bottom_hundredth / 30), top=bottom_hundredth)
                plt.xlim(0, 1)
                plt.rc('axes', linewidth=2)
                sns.despine(top=False, bottom=False, left=False, right=False)
                plt.savefig(out_dir + 'relative-NT-error-hist-zoom-' + str(sz).replace('.', '-') + '.' + ext,
                            dpi=300)
                plt.close('all')

                for r in runs:

                    tmp_dat = all_rel_pos.loc[all_rel_pos['Run'] == r]

                    if len(tmp_dat) == 0:
                        continue

                    fig = plt.figure(figsize=(sz, sz))
                    g = sns.displot(data=tmp_dat, x='Relative_Position', kind='hist', hue='Chain', height=sz,
                                    rug=True, rug_kws={"alpha": 0.3}, binwidth=bw, legend=False, hue_order=chains)
                    g.set_xlabels('Relative Position in TCR', fontweight='bold')
                    g.set_ylabels('Number Errors', fontweight='bold')
                    ax = g.axes[0, 0]
                    max_val = np.max(ax.get_ylim()[1])
                    plt.gca().set_ylim(bottom=-(max_val / 30))
                    plt.xlim(0,1)
                    plt.rc('axes', linewidth=2)
                    sns.despine(top=False, bottom=False, left=False, right=False)
                    plt.savefig(out_dir + 'relative-NT-error-hist_' + r.replace('\n(', '-').replace(')', '') + '-'
                                + str(sz).replace('.', '-') + '.' + ext, dpi=300)  #, bbox_inches='tight')
                    plt.close('all')

            if len(all_rel_pos_aa) > 0:
                used_runs = [x for x in cdr3_order_plotting if x in list(set(all_rel_pos_aa['Run']))]

                fig = plt.figure(figsize=(sz*4, sz))
                g = sns.displot(data=all_rel_pos_aa, x='Relative_Position', kind='hist', col='Run', hue='Chain',
                                height=sz, rug=True, rug_kws={"alpha": 0.3},
                                binwidth=bw, col_order=used_runs, hue_order=chains)
                g.set_xlabels('Relative Position in TCR', fontweight='bold')
                g.set_ylabels('Number Errors', fontweight='bold')
                axes = g.axes.flatten()
                for i in range(len(used_runs)):
                    axes[i].set_title(used_runs[i], fontweight='bold', verticalalignment='bottom')
                    axes[i].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
                    axes[i].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

                ax = g.axes[0, 0]
                max_val = np.max(ax.get_ylim()[1])
                plt.gca().set_ylim(bottom=-(max_val / 30), top=60)
                plt.xlim(0, 1)
                plt.rc('axes', linewidth=2)
                sns.despine(top=False, bottom=False, left=False, right=False)
                plt.savefig(out_dir + 'relative-AA-error-hist-' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
                plt.close('all')

            fxn.garbage_collection()

    plt.close('all')
    fxn.garbage_collection()

    return all_dat


def seamless_rename(in_df):
    """
    Change seamless run naming to better fit in a plot (without \n characters messing up the processing)
    :param in_df: Dataframe containing run names in the format "SLxxx" (where xxx = length of padding, 010/020/030/200)
    :return: Same dataframe with run names converted to "SL\n(xxx)" format
    """
    return in_df.\
        replace('SL200', 'SL\n(200)').\
        replace('SL010', 'SL\n(10)').\
        replace('SL020', 'SL\n(20)').\
        replace('SL030', 'SL\n(30)')


def tidy_dat_columns(df_to_tidy):
    """
    :param df_to_tidy: Pandas dataframe of TCR data
    :return: Same dataframe but with some additional columns/others renamed
    """
    df_to_tidy = df_to_tidy.replace('Benchmarking', 'Generated').replace('Processed', 'Original')
    df_to_tidy['Time (min)'] = df_to_tidy['Time (s)'] / 60
    df_to_tidy['% Stitched'] = (df_to_tidy['# Stitched'] / df_to_tidy['# TCRs']) * 100

    # Also produce the data specific names
    heather_dat = df_to_tidy.loc[df_to_tidy['Source'] == 'Heather']
    heather_dat['Provided CDR3'] = [x.replace('Heather_', '') for x in heather_dat['Sample']]

    sim_dat = df_to_tidy.loc[df_to_tidy['Source'] == 'immuneSIM']
    sim_dat['Provided CDR3'] = [x.replace('immuneSIM_', '') for x in sim_dat['Sample']]

    for x in ['10', '200', '20', '30']:
        heather_dat = heather_dat.replace('NT-SL-' + x, 'SL\n(' + x + ')')
        sim_dat = sim_dat.replace('NT-SL-' + x, 'SL\n(' + x + ')')

    return df_to_tidy, heather_dat, sim_dat


def plot_basic_stitch_stats(save_dir, plot_dat):
    """
    :param save_dir: Str path to directory to save basic plots in
    :param plot_dat: Pandas data frame containing information to be plotted
    :return: nothing
    """

    plt.rc('axes', linewidth=2)
    x_order = [seamless_rename(x) for x in cdr3_order]

    # Basic stats first
    for ext in fxn.exts:
        for sz in fxn.sizes:
            pad = 0.15

            # Timings
            fig = plt.figure(figsize=(sz, sz))
            g = sns.barplot(data=plot_dat, x='Provided CDR3', y='Time (min)', order=x_order, color='darkgray')
            for x in range(len(x_order)):
                row = plot_dat.loc[plot_dat['Provided CDR3'] == x_order[x]]
                y_val = row['Time (min)'].iloc[0]
                g.text(x, y_val + 0.1, round(row['Time (min)'].iloc[0], 1), color='black', ha="center", fontsize=14)

            sns.despine(top=True, bottom=False, left=False, right=True)
            plt.ylabel('Run time (min)', fontweight='bold')
            plt.xlabel('')
            plt.subplots_adjust(left=pad, bottom=pad)
            plt.savefig(save_dir + 'barplot-time-mins.' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            # and % stitched
            fig = plt.figure(figsize=(sz, sz))
            sns.barplot(data=plot_dat, x='Provided CDR3', y='% Stitched', order=x_order, color='darkgray')
            plt.ylim([99.96, 100])
            plt.axhline(y=100, alpha=.5, color='black', ls='--', lw=3)
            sns.despine(top=True, bottom=False, left=False, right=True)
            plt.ylabel('% TCRs stitched', fontweight='bold')
            plt.xlabel('')
            plt.subplots_adjust(left=pad, bottom=pad)
            plt.savefig(save_dir + 'barplot-pc-stitched.' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            fxn.garbage_collection()


# Set variables needed for this/replotting script
trac = 'ATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAA' \
       'GTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCAT' \
       'GTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGA' \
       'ACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGCTG'
trbc1 = 'AGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTC' \
        'TTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTC' \
        'CAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACG' \
        'AGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTG' \
        'TCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC'
trbc2 = 'AGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCT' \
        'ACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCA' \
        'GATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGT' \
        'GGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTCACCTCCGAGTCTTACCAGCAAGGGGTCCTGTCTG' \
        'CCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTCAAGAGAAAGGATTCCAGAGGC'

plt.rcParams.update({'figure.max_open_warning': 0})
base_out_dir = fxn.out_folder('plotting')
tcr_range = [(1*10)**exp for exp in range(2, 7)]
exp_tcr_range = [r"$10^{%1d}$" % x for x in range(2, 7)]
time_range_m = [(1*10)**x for x in [-3, -2, -1, 0, 1]]
cdr3_order = ['AA', 'NT', 'SL010', 'SL020', 'SL030', 'SL200']
cdr3_order_plotting = ['AA', 'NT', 'SL\n(20)', 'SL\n(200)']
cdr3_order_plotting_full = ['AA', 'NT', 'SL\n(10)', 'SL\n(20)', 'SL\n(30)', 'SL\n(200)']
error_bins = ['0', '1-2', '3-5', '6+']
res_order = ['NT', 'AA']
chains = ['TRA', 'TRB']
cdr3_palette = [sns.color_palette('bright')[x] for x in range(len(sns.color_palette('bright'))) if x in [4, 7, 9]]


if __name__ == "__main__":

    scripts_dir = fxn.check_scripts_cwd()

    plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial', 'font.weight': 'bold',
                         'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial',
                         'mathtext.bf': 'Arial:bold', 'mathtext.default': 'bf'})

    # Read in data, and rename some fields for more logical plotting names, and calculate some derived columns
    in_dat, stitched_dir = fxn.open_previous_data('stitched-results', 'timing_data', 'tsv')
    dat, heather, sim = tidy_dat_columns(in_dat)

    for chain in chains:
        dat.loc[(dat['Source'] == 'VDJdb') & (dat['Sample'].str.contains(chain)), 'Source'] = 'VDJdb ' + chain

    ####################################################################################################################

    # Then plot first bit of benchmarking: published sequences, Emerson AA/VDJdb AA

    # Pull out the published/third party data sources (VDJdb/Emerson)
    ord = ['Emerson', 'VDJdb', 'VDJdb TRA', 'VDJdb TRB']
    published = dat.loc[dat['Source'].isin(ord)]
    published_aa = published.loc[published['CDR3 input'] == 'AA']

    out_dir = fxn.sub_out_folder(base_out_dir, 'benchmarking-AA')

    for ext in fxn.exts:
        for sz in fxn.sizes:

            # Plot linear
            fig = plt.figure(figsize=(sz, sz))
            sns.lmplot(data=published_aa, x='# TCRs', y='Time (min)', hue='Source', markers=False, legend=False,
                       line_kws={'alpha': 0.5}, scatter_kws={'s': 0}, ci=None, truncate=True, lowess=True,
                       height=sz, aspect=1, hue_order=ord)
            sns.scatterplot(data=published_aa, x='# TCRs', y='Time (min)', hue='Source', style='Sampling',
                            alpha=.95, s=80, legend=False, markers=['o', '^'], hue_order=ord)
            plt.ylabel('Run time (min)', fontweight='bold')
            plt.xlabel('Number TCRs', fontweight='bold')
            ax = plt.gca()

            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            plt.xlim([90, 1.1e6])
            plt.yticks([0, 2, 4, 6, 8, 10])
            plt.xticks([0, 10, 100000, 500000, 1000000])
            plt.ticklabel_format(style='scientific', axis='x', useOffset=False)
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'correlation-linear-' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            # Plot log
            fig = plt.figure(figsize=(sz, sz))
            sns.lmplot(data=published_aa, x='# TCRs', y='Time (min)', hue='Source', markers=False, legend=False,
                       line_kws={'alpha': 0.5}, scatter_kws={'s': 0}, ci=None, truncate=True, lowess=True,
                       height=sz, aspect=1, hue_order=ord)
            sns.scatterplot(data=published_aa, x='# TCRs', y='Time (min)', hue='Source', style='Sampling',
                            alpha=.5, s=80, legend=False, markers=['o', '^'], hue_order=ord)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim([90, 1.2e6])
            plt.ylim([9e-4, 1.01e1])
            plt.xticks(tcr_range)
            plt.yticks(time_range_m)
            plt.ylabel('Run time (min)', fontweight='bold')
            plt.xlabel('Number TCRs', fontweight='bold')

            ax = plt.gca()
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            xlocmin = ticker.LogLocator(base=10.0, subs=[x / 10 for x in range(1, 10)], numticks=12)
            ax.xaxis.set_minor_locator(xlocmin)
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())

            ylocmin = ticker.LogLocator(base=10.0, subs=[x / 10 for x in range(1, 10)], numticks=12)
            ax.yaxis.set_minor_locator(ylocmin)
            ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            plt.rc('axes', linewidth=2)
            plt.subplots_adjust(left=0.15)
            plt.savefig(out_dir + 'correlation-log-' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            # And plot the legend
            fig = plt.figure(figsize=(sz, sz))
            g = sns.scatterplot(data=published_aa, x='# TCRs', y='Time (s)', hue='Source', style='Sampling',
                                alpha=0, s=80, markers=['o', '^'], hue_order=ord)
            g.set_xticklabels([''])
            g.set_yticklabels([''])
            g.set_xticks([])
            g.set_yticks([])
            plt.ylabel('')
            plt.xlabel('')
            sns.despine(top=True, bottom=True, left=True, right=True)
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'correlation-legend-' + str(sz).replace('.', '-') + '.' + ext,
                        dpi=300, bbox_inches='tight')
            plt.close('all')

            aa_mil = published_aa.loc[published_aa['# TCRs'] == 1e6]
            aa_mil = aa_mil.replace('VDJdb TRA', 'VDJdb\nTRA').replace('VDJdb TRB', 'VDJdb\nTRB')

            fig = plt.figure(figsize=(sz*1.35, sz))
            sns.barplot(data=aa_mil, x='Source', y='% Stitched')
            sns.swarmplot(data=aa_mil, x='Source', y='% Stitched', edgecolor='black', linewidth=1, size=7)
            plt.ylabel('% Stitched', fontweight='bold')
            plt.xlabel('')
            plt.ylim([99.9, 100.01])
            plt.subplots_adjust(left=0.24, bottom=0.17)
            sns.despine(top=True, bottom=False, left=False, right=True)
            plt.rc('axes', linewidth=2)
            plt.savefig(out_dir + 'percent-stitched-barplot-' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            fxn.garbage_collection()

    ####################################################################################################################

    # Then plot the next bit benchmarking: published sequences, Emerson AA/NT

    out_dir = fxn.sub_out_folder(base_out_dir, 'benchmarking-emerson-cdr3types')

    emerson = published.loc[(published['Source'] == 'Emerson') & (published['Sampling'] == 'Generated')]

    for ext in fxn.exts:
        for sz in fxn.sizes:
            fig = plt.figure(figsize=(sz*1.1, sz))
            sns.stripplot(data=emerson, x='# TCRs', y='Time (min)', hue='CDR3 input', order=tcr_range)

            ax = plt.gca()
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontweight('bold')

            plt.xticks(range(5), exp_tcr_range)
            # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.yscale('log')
            plt.ylabel('Run time (min)', fontweight='bold')
            plt.xlabel('Number TCRs', fontweight='bold')
            plt.rc('axes', linewidth=2)
            plt.subplots_adjust(left=0.23, bottom=0.18)
            sns.despine(top=False, bottom=False, left=False, right=False)
            plt.savefig(out_dir + 'stripplot.' + str(sz).replace('.', '-') + '.' + ext, dpi=300)
            plt.close('all')

            fxn.garbage_collection()

    ####################################################################################################################

    # Plot the basic benchmarking plots
    plot_basic_stitch_stats(fxn.sub_out_folder(base_out_dir, 'benchmarking-heather-cdr3types'), heather)
    plot_basic_stitch_stats(fxn.sub_out_folder(base_out_dir, 'benchmarking-immunesim-cdr3types'), sim)

    # Then the more elaborate ones
    heather_dat = plot_multi_modal(stitched_dir, fxn.proc_heather_dir, 'Heather', 'Heather',
                         'pre-stitchr-TCRs.tsv', 'inter_tag_seq')
    del heather_dat  # Deleting for memory conservation, comment as needed

    immunesim_dat = plot_multi_modal(stitched_dir, fxn.proc_immunesim_dir, 'immuneSIM', 'immuneSIM',
                         'pre-stitchr-TCRs.tsv', 'sequence')
    del immunesim_dat  # Deleting for memory conservation, comment as needed
