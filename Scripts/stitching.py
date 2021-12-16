# -*- coding: utf-8 -*-

"""
stitching.py

This script runs thimble (the high-throughput version of stitchr) on the files produced by preparation.py,
recording its timing as it does

"""

import os
import subprocess
import functions as fxn
from datetime import datetime

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.2.0'
__author__ = 'Jamie Heather'


def escape(path):
    """
    :param path: A string containing a path to be handled by bash
    :return: a suitably escaped version of that path
    """
    return path.replace(' ', '\\ ').replace('(', '\\(').replace(')', '\\)')


def run_thimble(list_of_dirs, output_dir_name, additional_genes, home_dir):
    """
    :param list_of_dirs: list of directories in which to look for files to stitch
    :param output_dir_name: name of directory that results will be put in
    :param additional_genes: boolean as to whether to use the -xg flag during stitching
    :param home_dir: home 'Scripts/' directory to return to
    :return:
    """

    abs_paths = []
    for d in list_of_dirs:
        print('\t' + d)
        in_files = [x for x in os.listdir(d) if '.tsv' in x and
                    ('bulk' in x or 'Heather' in x or 'HIP' in x or 'vdjdb' in x or 'immuneSIM' in x)]
        os.chdir(d)
        dir_path = os.getcwd() + '/'
        for f in in_files:
            abs_paths.append(dir_path + f)
        os.chdir(home_dir)

    abs_paths.sort()
    out_dir = fxn.out_folder(output_dir_name)
    os.chdir(out_dir)
    out_dir = os.getcwd() + '/'

    # And then move into the neighbouring stitchr directory to actually run the code
    os.chdir(home_dir + fxn.supp_script_dir + 'stitchr/Scripts/')

    dat = []
    print("Running Thimble...")
    for f in abs_paths:

        # Define the thimble command
        out_nam = f.split('/')[-1].split('.')[0]
        out_file = escape(out_dir + out_nam)
        print('\t' + out_nam)
        cmd = 'python3 thimble.py -in ' + escape(f) + ' -o ' + out_file

        # Add additional flags
        if '-SL' in out_nam:
            cmd += ' -sl'

        if additional_genes:
            cmd += ' -xg'

        # Start a timer, then run that command
        start_time = datetime.now()
        subprocess.call(cmd, shell=True)
        time_taken = datetime.now() - start_time
        bits = f.split('/')

        # Define some of the fields of the output data table
        if 'SL' in out_nam:
            cdr = 'SL'
        elif 'AA' in out_nam:
            cdr = 'AA'
        else:
            cdr = 'NT'

        if 'bulk' in out_nam:
            bulk_bits = out_nam.split('_')
            rep = bulk_bits[5]
        else:
            rep = ''

        tcrs = 0  # Total number of TCRs input...
        stitched = 0  # ... and how many of those produced a stitched TCR
        out_path = out_dir + out_nam + '.tsv'
        with fxn.opener(out_path, 'r') as in_file:
            line_count = 0
            for line in in_file:
                if line_count > 0:  # Don't count header row
                    tcrs += 1
                    line_bits = line.split('\t')
                    if fxn.dna_check(line_bits[1]) or fxn.dna_check(line_bits[2]):
                        stitched += 1
                line_count += 1

        # Output not yet zipped (to get accurate timing), so run that now (for space/downstream file expectations)
        subprocess.check_call(['gzip', out_path])  # TODO add force flag (-f)?

        fxn.garbage_collection()  # Force garbage collection to try to increase fairness of timing between iterations

        dat.append([bits[-3], bits[-2], cdr, out_nam, tcrs, stitched, rep, time_taken.total_seconds()])

    headers = ['Source', 'Sampling', 'CDR3 input', 'Sample', '# TCRs', '# Stitched', 'Repeat', 'Time (s)']

    dat = fxn.list_to_df(dat, headers, False)

    dat.to_csv(out_dir + 'timing_data.tsv', sep='\t', index=False)

    os.chdir(home_dir)
    return dat

    # TODO check gzip works (and deleted old one)


# Sort out the directory business:
# define absolute paths to all of theoutput_dir_nameoutput_dir_name input files and the the output folder
if __name__ == "__main__":
    scripts_dir = fxn.check_scripts_cwd()

    print("Finding input data...")
    in_dirs = [fxn.proc_heather_dir, fxn.proc_emerson_dir, fxn.proc_vdjdb_dir, fxn.proc_immunesim_dir,
               fxn.bench_emerson_dir, fxn.bench_vdjdb_dir]

    thimble_dat = run_thimble(in_dirs, 'stitched-results', False, scripts_dir)
