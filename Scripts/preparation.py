# -*- coding: utf-8 -*-

"""
preparation.py

Takes the raw data provided in this repo, and automates production of the input files required for downstream analysis

"""

import os
import shutil
import collections as coll
import functions as fxn
import pandas as pd
import random as rnd
from copy import deepcopy
from functions import *

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.0'
__author__ = 'Jamie Heather'


def get_padded_junction(wider_seq, seq_to_pad, pad_len):
    """
    Fxn to take a given long read TCR plus it's CDR3 junction sequence, and provide that CDR3 plus a certain context
    :param wider_seq: the wider sequence context from which context will be drawn (i.e. the inter tag sequence)
    :param seq_to_pad: the substring sequence to pad (i.e. the CDR3 junction)
    :param pad_len: integer number of nucleotides to garner from either side of the seq_to_pad from the wider_seq
    :return: sequence corresponding to substring extended by the appropriate number of characters from the wider str
    """
    substr_index = wider_seq.index(seq_to_pad)
    start_pos = substr_index - pad_len
    if start_pos < 0:
        start_pos = 0
    end_pos = substr_index + pad_len + len(seq_to_pad)
    return wider_seq[start_pos:end_pos]


def sample_repertoires(list_of_tcrs, out_dir, sample_prefix, file_headers):
    """
    :param list_of_tcrs: list of lines read in from an as-is thimble input file
    :param out_dir: path to directory to save results in
    :param sample_prefix: name of sample to include in the individual file names
    :param file_headers: first row to write to the outfile
    :return: nothing, but the tcr list gets sub/up-sampled across a range of sizes and size
    """
    for exp in range(2, 7):
        size = (1*10)**exp
        for repeat in range(1, 4):
            out_nam = out_dir + sample_prefix + \
                      str(size) + '_rep_' + str(repeat) + '.tsv'
            with open(out_nam, 'w') as out_file:
                tcrs = rnd.choices(list_of_tcrs, k=size)
                out_file.write(file_headers + ''.join(tcrs))


def airr_to_stitchr(in_df, out_dir, sample_nam, stitchr_file_headers, long_seq_field):
    """
    Convert AIRR-C standard TCR files into appropriate Thimble-compatible files
    :param in_df: Pandas dataframe of AIRR-C standard data
    :param out_dir: Str path to save out data
    :param sample_nam: Str name for sample
    :param stitchr_file_headers: List of strings relating to output columns for Thimble input
    :param long_seq_field: Str name of field containing the full length TCR sequence to be used
    :return:
    """

    # Loop through the retained data, slotting into the appropriate place in the stitchr format table
    # Output matrices of tables for thimble input, giving the same TCRs with...:
    out_aa = []  # amino acid CDR3 junctions
    out_nt = []  # nucleotide CDR3 junctions
    out_nt_ext = coll.defaultdict(list)

    conversion = {
        'V': 'v_call',
        'J': 'j_call',
        '_CDR3': 'junction'
    }

    for row in in_df.index:
        tmp_dat = in_df.loc[row]
        conv_row = [''] * len(stitchr_file_headers)
        conv_row[0] = str(row)  # name the TCR based on row number, make tracking TCRs easier downstream

        # Pull out the chain from the J gene
        chain = tmp_dat['j_call'][:3]
        for field in conversion:
            out_id = chain + field
            idx = stitchr_file_headers.index(out_id)
            conv_row[idx] = tmp_dat[conversion[field]]

        # And also add in a simple ATG leader (as stitchr requires something)
        conv_row[0] = int(conv_row[0])
        conv_row[stitchr_file_headers.index(chain + '_leader')] = 'ATG'

        # Output the simple already available nt CDR3 data
        out_nt.append(deepcopy(conv_row))

        # Then generate each of the others and store those
        conv_row[idx] = fxn.translate_nt(tmp_dat['junction'])
        out_aa.append(deepcopy(conv_row))
        for i in [10, 20, 30, 200]:
            conv_row[idx] = get_padded_junction(tmp_dat[long_seq_field], tmp_dat['junction'], i)
            out_nt_ext[i].append(deepcopy(conv_row))

    out_aa = fxn.list_to_df(out_aa, stitchr_file_headers, False)
    out_nt = fxn.list_to_df(out_nt, stitchr_file_headers, False)

    out_aa.to_csv(out_dir + sample_nam + '_AA.tsv.gz', compression='gzip', sep='\t', index=False)
    out_nt.to_csv(out_dir + sample_nam + '_NT.tsv.gz', compression='gzip', sep='\t', index=False)

    for i in out_nt_ext:
        out_nt_ext[i] = fxn.list_to_df(out_nt_ext[i], stitchr_file_headers, False)
        out_nt_ext[i].to_csv(out_dir + sample_nam + '_NT-SL-' + str(i) + '.tsv.gz', compression='gzip',
                             sep='\t', index=False)


def format_heather_data(path_to_file):
    """
    :param path_to_file: Path to autoDCR output file of rearranged TCRs from the UMI-corrected Heather et al data
    :return: pandas dataframe containing filtered TCRs with only necessary info, in the right format
    """

    dat = pd.read_csv(path_to_file, sep='\t', compression='infer')
    dat = dat.dropna(axis=1, how='all')
    rearranged = dat.loc[dat['productive'] == 'T']

    # Then collapse redundant rows, keeping only the info we care about
    collapsed = rearranged['sequence_id'].groupby(
        [rearranged.v_call, rearranged.j_call,
         rearranged.junction, rearranged.inter_tag_seq,
         rearranged.v_jump, rearranged.j_jump]
        ).apply(list).reset_index()

    # Remove any entries that lack unambiguous allelic determination, or have very short/very long CDR3s
    non_ambig = collapsed.loc[~(collapsed['v_call'].str.contains(',')) &
                              ~(collapsed['j_call'].str.contains(',')) &
                              ~(collapsed['junction'].str.len() < 8 * 3) &
                              ~(collapsed['junction'].str.len() > 25 * 3)]

    # Let's count the number of donors each occurs in
    non_ambig['number_donors'] = [len(set([y[:4] for y in x])) for x in non_ambig['sequence_id']]

    # And let's filter out the tiny number of trans-locus (i.e. TRAV-TRBJ/TRBV-TRAJ) probable PCR chimeras
    non_ambig = non_ambig.loc[((non_ambig['v_call'].str.contains('TRA')) & (non_ambig['j_call'].str.contains('TRA'))) |
                              ((non_ambig['v_call'].str.contains('TRD')) & (non_ambig['j_call'].str.contains('TRA'))) |
                              ((non_ambig['v_call'].str.contains('TRB')) & (non_ambig['j_call'].str.contains('TRB'))) ]

    # Also filter out the tiny minority of rearrangements (1 or 2) which don't have the junction detectable in the ITS
    non_ambig = non_ambig[non_ambig.apply(lambda x: x.junction in x.inter_tag_seq, axis=1)]

    # And then let's focus on those where we have full V genes,
    # i.e. the v_jump (point in the full V gene the most upstream tag was found) == 0
    full_tcrs = non_ambig.loc[non_ambig['v_jump'] == 0]

    # Also need to get rid of any TCRs found to be using genes for which IMGT only has partial information
    # (as further edits would be required to align the inter-tag sequence with the stitched!)
    partial_genes = []
    with open('Supplementary-Scripts/autoDCR/human.fasta', 'r') as germlines:
        for line in germlines:
            if 'partial' in line:
                partial_genes.append(line.split('|')[1])

    full_tcrs = full_tcrs.loc[~full_tcrs['v_call'].isin(partial_genes)]
    return full_tcrs



if __name__ == "__main__":

    fxn.check_scripts_cwd()

    # First, check that all necessary directories are present, and create them if they're not
    all_dirs = [x for x in dir() if x.endswith('_dir')]
    for d in all_dirs:
        fxn.check_directory(vars()[d])

    # Download all supplementary scripts via git
    os.chdir(fxn.supp_script_dir)
    for git in fxn.gits:
        # Start with a fresh directory, just in case
        if git in os.listdir(os.getcwd()):
            fxn.run_bash('rm -rf ' + git)

        fxn.run_bash('git clone ' + fxn.gits[git])
        # TODO could here add a command that cds into the script dir and checks out a branch with right versions

    os.chdir('../')

    # Emerson data
    print("Processing Emerson (Adaptive) data...")
    print("\tConverting to AIRR-C standard")

    #  First convert to AIRR Community file formats/IMGT gene nomenclature
    raw_emerson_files = [x for x in os.listdir(fxn.raw_emerson_dir) if x.startswith('HIP') and '.tsv' in x]

    for raw in raw_emerson_files:
        nam = raw.split('.')[0]
        fixed_nam = nam + '-fixed.tsv.gz'
        print("\t\t" + nam)

        # We actually need to first correct for a feature of the Emerson data,
        # in which allele information is stored in its own different column = copy across to the V/J gene columns
        with fxn.opener(fxn.raw_emerson_dir + raw, 'r') as in_file, \
                fxn.opener(fxn.int_emerson_dir + fixed_nam, 'w') as out_file:
            line_count = 0
            for line in in_file:
                save = True
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    headers = bits
                else:
                    for gene in ['v', 'j']:
                        gene_call = bits[headers.index(gene + '_gene')]
                        allele_call = bits[headers.index(gene + '_allele')]

                        # Filter out ambiguous gene entries
                        if ',' in gene_call:
                            save = False
                            continue

                        # ... and if there's an allele call, add that information to the gene level column
                        elif allele_call:
                            if not allele_call.startswith('0'):
                                allele_call = '0' + allele_call
                            bits[headers.index(gene + '_gene')] = gene_call + '*' + allele_call

                if save:
                    out_file.write('\t'.join(bits) + '\n')

                line_count += 1

        # Ignore: D information / orphon genes / non-productive reads / non C-F junctions / CDR3s < 8 residues
        fxn.run_bash("python3 Supplementary-Scripts/immunoseq2airr/immunoseq2airr.py -nd -a -or -pf -mf -z -jlf 8 " \
              "-i " + fxn.int_emerson_dir + fixed_nam + " " \
              "-o " + fxn.int_emerson_dir + nam + "-airr " \
              "-p " + "Supplementary-Scripts/immunoseq2airr/emerson-parameters.tsv")

        os.remove(fxn.int_emerson_dir + fixed_nam)

    # Then process those into thimble files with either amino acid or entire read nt CDR3 sequences
    airr_emerson_files = [x for x in os.listdir(fxn.int_emerson_dir) if x.startswith('HIP') and '-airr.tsv' in x]

    emerson_conversion = {
        'TRBV': 'v_call',
        'TRBJ': 'j_call',
        'TRB_CDR3': 'junction_aa',
        'TCR_name': 'sequence_id'
    }

    print("\tConverting to Thimble format")
    for airr in airr_emerson_files:
        nam = airr.split('-')[0]
        print("\t\t" + nam)
        df = pd.read_csv(fxn.int_emerson_dir + airr, sep='\t')
        df = df.loc[(~df['v_call'].str.contains(',')) & (~df['j_call'].str.contains(','))]

        # Open outfiles, add headers
        with fxn.opener(fxn.proc_emerson_dir + nam + '_AA.tsv.gz', 'w') as aa_out, \
                fxn.opener(fxn.proc_emerson_dir + nam + '_NT-SL.tsv.gz', 'w') as sl_out:
            aa_out.write('\t'.join(fxn.stitchr_headers) + '\n')
            sl_out.write('\t'.join(fxn.stitchr_headers) + '\n')

            # Then go through the filtered dataframe, writing out the TCRs with either the amino acid or whole nt CDR3
            for tcr in df.index:
                out_row = [''] * len(fxn.stitchr_headers)
                for field in emerson_conversion:
                    out_row[fxn.stitchr_headers.index(field)] = df.loc[tcr][emerson_conversion[field]]

                aa_out.write('\t'.join(out_row) + '\n')
                out_row[fxn.stitchr_headers.index('TRB_CDR3')] = df.loc[tcr]['sequence']

                sl_out.write('\t'.join(out_row) + '\n')

    print("\tRandomly sampling for benchmarking")
    # We need to separately pool the AA and NT-SL results across all 5 files, then randomly sample from those

    aa_emerson_files = [x for x in os.listdir(fxn.proc_emerson_dir) if x.startswith('HIP') and '_AA.tsv' in x]
    sl_emerson_files = [x for x in os.listdir(fxn.proc_emerson_dir) if x.startswith('HIP') and '_NT-SL.tsv' in x]

    benchmarking_dat = coll.defaultdict(list)

    conv = {'aa_': 'AA_', 'sl_': 'NT-SL_'}

    # First read the data in to lists
    for input_type in ['aa_', 'sl_']:
        for f in vars()[input_type + 'emerson_files']:
            with fxn.opener(fxn.proc_emerson_dir + f, 'r') as in_file:
                line_count = 0
                for line in in_file:
                    if line_count == 0:
                        headers = line
                    else:
                        benchmarking_dat[input_type + 'emerson_in'].append(line)
                    line_count += 1
        # Then randomly up/downsample those to particular file sizes
        sample_repertoires(benchmarking_dat[input_type + 'emerson_in'], fxn.bench_emerson_dir,
                           'bulk_emerson_' + conv[input_type], headers)

    # Then VDJdb data
    print("Processing VDJdb data...")
    print("\tConverting to Thimble format")
    vdjdb_conversion = {
        'V': 7,
        'J': 8,
        '_CDR3': 1
    }

    # Write out a file of all TCRs, plus separate alpha/beta chains
    with fxn.opener(fxn.raw_vdjdb_dir + 'vdjdb.slim.txt', 'r') as in_file, \
            fxn.opener(fxn.proc_vdjdb_dir + 'vdjdb_AA_all.tsv.gz', 'w') as out_file, \
            fxn.opener(fxn.proc_vdjdb_dir + 'vdjdb_AA_TRA.tsv.gz', 'w') as tra_out_file, \
            fxn.opener(fxn.proc_vdjdb_dir + 'vdjdb_AA_TRB.tsv.gz', 'w') as trb_out_file:
        out_file.write('\t'.join(fxn.stitchr_headers) + '\n')
        tra_out_file.write('\t'.join(fxn.stitchr_headers) + '\n')
        trb_out_file.write('\t'.join(fxn.stitchr_headers) + '\n')
        line_count = 0
        for line in in_file:
            bits = line.rstrip().split('\t')
            if line_count == 0:
                headers = bits

            # Select human TCRs with V/J/CDR3 information
            elif bits[2] == 'HomoSapiens' and bits[7] and bits[8] and bits[1]:

                # Again, get rid of ambiguous gene calls
                if ',' in bits[7] or ',' in bits[8]:
                    continue

                # Screen out dubious CDR3s
                if len(bits[1]) < 8 or bits[1][-1] not in ['F', 'W']:
                    continue

                # Then output to the appropriate columns
                chain = bits[0]
                out_row = [''] * len(fxn.stitchr_headers)
                out_row[0] = 'VDJdb|' + str(line_count).zfill(5)

                for field in ['V', 'J', '_CDR3']:
                    out_row[fxn.stitchr_headers.index(chain + field)] = bits[vdjdb_conversion[field]]

                out_file.write('\t'.join(out_row) + '\n')

                if chain == 'TRA':
                    tra_out_file.write('\t'.join(out_row) + '\n')
                elif chain == 'TRB':
                    trb_out_file.write('\t'.join(out_row) + '\n')
                else:
                    raise ValueError("Rogue chain detected!")

            line_count += 1

    print("\tRandomly sampling for benchmarking")
    # And then again, duplicate up and sample back to fixed sizes

    # First read the data in to lists
    for f in [x for x in os.listdir(fxn.proc_vdjdb_dir) if x.startswith('vdjdb_AA') and '.tsv' in x]:
        out_nam = f.replace('tsv.', '').replace('.gz', '').replace('_AA', '').replace('_', '-')
        with fxn.opener(fxn.proc_vdjdb_dir + f, 'r') as in_file:
            line_count = 0
            for line in in_file:
                if line_count == 0:
                    headers = line
                else:
                    benchmarking_dat[out_nam].append(line)
                line_count += 1
        # Then randomly up/downsample those to particular file sizes
        sample_repertoires(benchmarking_dat[out_nam], fxn.bench_vdjdb_dir,
                           'bulk_' + out_nam + '_AA_', headers)

    del benchmarking_dat

    # And finally the UMI-barcoded merged reads from Heather et al
    print("Processing Heather et al data...")
    print("\tDownloading merged FASTQ from SRA")

    # First download the raw data from SRA (using the entries as stored in the README table)
    with fxn.opener('../README.md', 'r') as in_file:
        for line in in_file:
            if 'SRR' in line:
                bits = line.rstrip().split(' ')
                donor = bits[1]
                chain = bits[3]
                srr = bits[7]
                out_nam = donor + '_' + chain + '.fastq.gz'
                print('\t\t\t' + srr + '\t|\t' + out_nam)
                fxn.run_bash("prefetch " + srr)
                fxn.run_bash("fastq-dump " + srr + " --gzip --readids --outdir " + fxn.raw_heather_dir)

                # Rename the downloaded FASTQ, delete the pre-fetched SRA
                os.rename(fxn.raw_heather_dir + srr + '.fastq.gz', fxn.raw_heather_dir + out_nam)
                shutil.rmtree(srr)

    # Then go through these FASTQs and stringently collapse into high-quality full-length read FASTA files

    # Find all fastq files (that meet my format)
    all_fq = [x for x in os.listdir(fxn.raw_heather_dir) if x.startswith('HV')
              and (x.endswith(".fastq.gz") or x.endswith(".fastq"))]
    all_fq.sort()

    print("\tCollapsing merged FASTQ into high-quality FASTA")
    if len(all_fq) == 0:
        raise IOError("No FASTQ files (*fq/*fq.gz) detected; nothing to collapse.")

    for fl in all_fq:
        print("\t\t" + fl)
        # Collapse all reads, returning a counter
        col_reads, counts = fxn.collapse(fxn.raw_heather_dir + fl)

        if fl.endswith(".fastq"):
            base_name = fl[:-6]
        elif fl.endswith(".fastq.gz"):
            base_name = fl[:-9]
        else:
            raise IOError("Undetermined FASTQ suffix information: " + fl)

        out_name = base_name + '.fasta.gz'

        # Write these collapsed transcripts out into their own fasta files
        fxn.write_fasta(col_reads, base_name.replace('_', '-'), fxn.int_heather_dir + out_name)

    # Then merge all these into one big fasta file
    fxn.run_bash('gunzip -c ' + fxn.int_heather_dir + 'HV*fasta.gz | gzip > ' + \
          fxn.int_heather_dir + 'HV_merged_combined.fasta.gz')

    # And assign TCRs with autoDCR
    print("\tAnnotating rearranged TCRs using autoDCR...")
    fxn.run_bash("python3 Supplementary-Scripts/autoDCR/autoDCR.py -fq " + fxn.int_heather_dir \
          + "HV_merged_combined.fasta.gz -o " + fxn.int_heather_dir + " -dd Supplementary-Scripts/ -or forward")

    print("\tConverting to Thimble format")
    full = format_heather_data(fxn.int_heather_dir + 'HV_merged_combined.tsv.gz')
    full.to_csv(fxn.proc_heather_dir + 'pre-stitchr-TCRs.tsv.gz', compression='gzip', sep='\t')
    airr_to_stitchr(full, fxn.proc_heather_dir, 'Heather', fxn.stitchr_headers, 'inter_tag_seq')

    # Also combine and convert the immuneSIM simulated TCR data
    raw_immunesim_files = [x for x in os.listdir(fxn.raw_immunesim_dir) if '.tsv' in x and '~' not in x]

    tmp_is_dfs = []
    for f in raw_immunesim_files:
        tmp_is_dfs.append(pd.read_csv(fxn.raw_immunesim_dir + f, compression='infer', sep='\t'))

    is_dat = pd.concat(tmp_is_dfs, ignore_index=True)

    # Throw out simulated TCRs which don't end in feasibly functional motifs
    is_dat = is_dat.loc[((is_dat['junction_aa'].str.endswith('F')) |
                         (is_dat['junction_aa'].str.endswith('W')) |
                         (is_dat['junction_aa'].str.endswith('C')))]

    is_dat.to_csv(fxn.proc_immunesim_dir + 'pre-stitchr-TCRs.tsv.gz', compression='gzip', sep='\t')

    airr_to_stitchr(is_dat, fxn.proc_immunesim_dir, 'immuneSIM', fxn.stitchr_headers, 'sequence')
