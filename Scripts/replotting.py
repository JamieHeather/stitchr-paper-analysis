# -*- coding: utf-8 -*-

"""
replotting.py

Re-perform basic plotting after inferring novel allele in the Heather dataset

"""

import os
import collections as coll
import functions as fxn
import preparation as prep
import plotting as plot
import stitching as stitch
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functions import *

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.2'
__author__ = 'Jamie Heather'

db_order = ['IMGT', 'IMGT+inferred']


if __name__ == "__main__":

    scripts_dir = fxn.check_scripts_cwd()

    # Check re-run directory exists
    fxn.check_directory(fxn.rerun_heather_dir)

    # Rerun autoDCR with the modified reference TCR files
    print("\tAnnotating rearranged TCRs using autoDCR...")
    fxn.run_bash("python3 " + fxn.supp_script_dir + "autoDCR/autoDCR.py -fq " + fxn.int_heather_dir +
                 "HV_merged_combined.fasta.gz -o " + fxn.rerun_heather_dir +
                 " -dd " + fxn.supp_script_dir + "autoDCR/ -jv -or forward -sp human-plus")

    print("\tConverting to Thimble format")

    full = prep.format_heather_data(fxn.rerun_heather_dir + 'HV_merged_combined.tsv.gz')
    full.to_csv(fxn.rerun_heather_dir + 'pre-stitchr-TCRs.tsv.gz', compression='gzip', sep='\t')
    prep.airr_to_stitchr(full, fxn.rerun_heather_dir, 'Heather', fxn.stitchr_headers, 'inter_tag_seq')

    print("\tRunning thimble on re-processed data")
    time_dat = stitch.run_thimble([fxn.rerun_heather_dir], 'stitched-rerun-results', True, scripts_dir)

    plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial', 'font.weight': 'bold',
                         'mathtext.fontset': 'custom', 'mathtext.it': 'Arial:italic', 'mathtext.rm': 'Arial',
                         'mathtext.bf': 'Arial:bold', 'mathtext.default': 'bf'})

    db_palette = [(0.30, 0.52, 0.66), (0.75, 0.15, 0.031)]

    in_dat, stitched_dir = fxn.open_previous_data('stitched-rerun-results', 'timing_data', 'tsv')
    tidied_dat, renamed_dat, unused_sim_dat = plot.tidy_dat_columns(in_dat)
    del tidied_dat
    del unused_sim_dat
    fxn.garbage_collection()

    plot.plot_basic_stitch_stats(fxn.sub_out_folder(plot.base_out_dir, 'benchmarking-heather-REPROCESS'), renamed_dat)

    rerun_dat = plot.plot_multi_modal(stitched_dir, fxn.rerun_heather_dir, 'Heather', 'Heather-HumanPlusReference',
                                      'pre-stitchr-TCRs.tsv', 'inter_tag_seq')

    # Then read in the most recent matching TCR data that was run with the baseline (non-inferred) TCR V/J genes
    last_dat = [x for x in os.listdir(fxn.output_dir)
                if x.endswith('-plotting')
                and 'Heather' in os.listdir(fxn.output_dir + x)]
    last_dat.sort()

    dat_file = 'all_dat.tsv.gz'
    for x in last_dat[::-1]:
        dir_path = fxn.output_dir + x + '/Heather/'
        dat_path = dir_path + dat_file
        print(x)
        if dat_file in os.listdir(dir_path):
            print('\tFound data file.')
            dat = pd.read_csv(dat_path, sep='\t', compression='infer', index_col=0)
            break
        else:
            print('\t\tFailed to find data file.')

    # Pull out TCRs that use an inferred potentially novel allele
    inferred_allele_tcr_dat = rerun_dat.loc[rerun_dat['V'].str.contains('_')]
    i_a_tcr_ids = [int(x) for x in set(inferred_allele_tcr_dat['TCR_ID'])]
    i_a_tcr_seqs = list(set(inferred_allele_tcr_dat['Sequence']))

    runs = ['AA', 'NT', 'SL\n(10)', 'SL\n(20)', 'SL\n(30)', 'SL\n(200)']

    # Find the TCRs shared with the previous analysis
    matched_tcrs_orig = plot.seamless_rename(dat.loc[dat['Sequence'].isin(i_a_tcr_seqs)])
    matched_tcrs_rerun = rerun_dat.loc[rerun_dat['Sequence'].isin(list(set(matched_tcrs_orig['Sequence'])))]
    matched_tcrs = list(set(matched_tcrs_orig['Sequence']))

    # TCR_IDs aren't maintained successfully across multiple analyses
    # We therefore need to go through sequence-by-sequence (which are unique)
    comparison_dat = []
    convert_type = {'Edit_Distance': 'NT', 'Edit_Distance_AA': 'AA'}
    for tcr in matched_tcrs:
        for run in runs:
            out_line_base = [tcr, run]
            tmp_orig_dat = matched_tcrs_orig.loc[(matched_tcrs_orig['Sequence'] == tcr) &
                                                 (matched_tcrs_orig['Run'] == run)].iloc[0]
            tmp_rerun_dat = matched_tcrs_rerun.loc[(matched_tcrs_rerun['Sequence'] == tcr) &
                                                   (matched_tcrs_rerun['Run'] == run)].iloc[0]
            for edit in convert_type:
                comparison_dat.append(out_line_base + [tmp_orig_dat[edit], convert_type[edit], 'IMGT',
                                      tmp_orig_dat['Number_' + convert_type[edit] + '_Correct'] / tmp_orig_dat[
                                          'Number_' + convert_type[edit]] * 100])
                comparison_dat.append(out_line_base + [tmp_rerun_dat[edit], convert_type[edit], 'IMGT+inferred',
                                      tmp_rerun_dat['Number_' + convert_type[edit] + '_Correct'] / tmp_rerun_dat[
                                          'Number_' + convert_type[edit]] * 100])

    cdat_headers = ['TCR', 'Run', 'Edit Distance', 'Residue', 'TCR database', '% Residues Correct']
    comparison_dat = fxn.list_to_df(comparison_dat, cdat_headers, False)

    out_dir = fxn.sub_out_folder(plot.base_out_dir, 'inferred-TCR-comparison')

    g = sns.catplot(data=plot.seamless_rename(comparison_dat), x='Run', y='Edit Distance', kind='bar',
                    col='Residue', hue='TCR database', palette=db_palette,
                    order=runs, hue_order=db_order, col_order=plot.res_order)

    count = 0
    for axes in g.axes.flat:
        axes.set_xlabel('')
        axes.set_title(axes.get_title(), fontweight='bold')
        if count == 0:
            axes.set_ylabel('Edit distance', fontweight='bold')
        count += 1

    plt.rc('axes', linewidth=2)
    plt.savefig(out_dir + 'edit-distances.png', dpi=300, bbox_inches='tight')
    plt.close()

    g = sns.catplot(data=plot.seamless_rename(comparison_dat), x='Run', y='% Residues Correct', kind='bar',
                    col='Residue', hue='TCR database', palette=db_palette,
                    order=runs, hue_order=db_order, col_order=plot.res_order)

    count = 0
    for axes in g.axes.flat:
        axes.set_xlabel('')
        axes.axhline(y=100, alpha=.8, color=plot.cdr3_palette[count], ls='--', lw=5)
        axes.set_title(axes.get_title(), fontweight='bold')
        if count == 0:
            axes.set_ylabel('Identical residues (%)', fontweight='bold')
        count += 1

    plt.ylim([98.5, 100])
    plt.rc('axes', linewidth=2)
    plt.savefig(out_dir + 'pc-residues-correct.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Calculate percent perfect
    perfect_dat = []
    for run in runs:
        for res in plot.res_order:
            for tcrdb in db_order:
                tmp_dat = comparison_dat.loc[(comparison_dat['Residue'] == res) &
                                             (comparison_dat['TCR database'] == tcrdb) &
                                             (comparison_dat['Run'] == run)]
                perfect_dat.append([run, res, tcrdb, len(tmp_dat.loc[tmp_dat['Edit Distance'] == 0])/len(tmp_dat)*100])

    perfect_dat = fxn.list_to_df(perfect_dat, ['Run', 'Residue', 'TCR database', '% Perfect'], False)
    perfect_dat = plot.seamless_rename(perfect_dat)

    offsets = [-0.185, 0.215]
    for sz in fxn.sizes:
        fig = plt.figure(figsize=(sz*2, sz))
        g = sns.catplot(data=perfect_dat, x='Run', y='% Perfect', col='Residue', hue='TCR database',
                        order=runs, kind='bar', hue_order=db_order, col_order=plot.res_order,
                        palette=db_palette, height=sz)

        count = 0
        for ax in g.axes.flat:

            # Only the leftmost plot needs a y axis label
            ax.set_xlabel('')
            if count == 0:
                ax.set_ylabel('Perfectly matching TCRs (%)', fontweight='bold')
            ax.set_title(ax.get_title().split(' ')[-1], fontweight='bold', verticalalignment='bottom')
            count += 1

            # Add text numbers above bars to make it easier to read off values
            res = ax.get_title().split(' ')[-1]
            for x in range(len(runs)):
                row = perfect_dat.loc[perfect_dat['Run'] == runs[x]]
                for db_pos in range(len(db_order)):
                    db = db_order[db_pos]
                    y_val = row.loc[(row['Residue'] == res) & (row['TCR database'] == db)]['% Perfect'][0]
                    if y_val < 20:
                        plot_y = y_val + 5
                        txt_col = db_palette[db_pos]
                    else:
                        plot_y = y_val - 15
                        txt_col = 'white'
                    ax.text(x + offsets[db_pos], plot_y, round(y_val, 1), color=txt_col,
                            ha="center", fontsize=14, rotation='vertical', horizontalalignment='center')

        plt.ylim([0, 100])
        plt.rc('axes', linewidth=2)
        plt.savefig(out_dir + 'pc-tcrs-perfect-' + str(sz).replace('.', '-') + '.png', dpi=300)#, bbox_inches='tight')
        plt.close()

    # Find the sequences that are keeping us shy of 100% accuracy!
    persistent_errors = comparison_dat.loc[
        (comparison_dat['Run'] == 'SL200') &
        (comparison_dat['Residue'] == 'AA') &
        (comparison_dat['TCR database'] == 'IMGT+inferred') &
        (comparison_dat['% Residues Correct'] != 100)]
