# -*- coding: utf-8 -*-

"""
functions.py

Provides shared functions and other objects for running the scripts in this repo,
which in turn exists to run and plot analyses testing out the stitchr suite of TCR sequence generation tools

"""

import datetime
import gzip
import os
import subprocess
import collections as coll
import Levenshtein as lev
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.0'
__author__ = 'Jamie Heather'


def check_scripts_cwd():
    """
    Check we're in the right directory (Scripts)
    :return: nothin
    """
    if not os.getcwd().endswith('/Scripts'):
        if 'Scripts' in os.listdir(os.getcwd()):
            os.chdir('Scripts')
        else:
            raise IOError("Check your working directory: this is designed to be run from the parent or Scripts folders")


def opener(file_path, open_mode):
    """
    :rtype: object
    :param file_path: path to file to be opened
    :param open_mode: mode by which to open the file (e.g. w/r/a)
    :return: the appropriate file opening command (open or gzip.open)
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, open_mode + 't')
    else:
        return open(file_path, open_mode)


def nest():
    """
    Create nested defaultdicts
    """
    return coll.defaultdict(list)


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def list_to_df(input_list, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param input_list: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(input_list)
    df = df.rename(index=str, columns=dict(zip(range(len(headers)), headers)))
    df = df.sort_values(by=[headers[0]])
    if rename is True:
        df = df.set_index(headers[0], drop=True)
    return df


def check_directory(dir_path):
    """
    Given the path to a directory, checks whether it exists: if it doesn't, it creates it
    :param dir_path: Path to be checked/created
    :return: Path to confirmed/created directory
    """
    if not dir_path.endswith('/'):
        dir_path = dir_path + '/'
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path


def translate_nt(nt_seq):
    """
    :param nt_seq: Nucleotide sequence to be translated
    :return: corresponding amino acid sequence
    """
    aa_seq = ''
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3].upper()
        if len(codon) == 3:
            try:
                aa_seq += codons[codon]
            except:
                raise IOError("Cannot translate codon: " + codon)
    return aa_seq


def rev_comp(seq):
    """
    :param seq: Any DNA string, composed of just the four bases (upper or lower case)
    :return: The reverse complement of that DNA string (maintaining same case)
    """
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def readfq(fp):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    :param fp: opened file containing fastq or fasta reads
    :yield: read id, read sequence, and (where available) read quality scores
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:], [], None  # This version takes the whole line (post '>')
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break

        else:  # this is a fastq record
            sequence, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(sequence):  # have read enough quality
                    last = None
                    yield name, sequence, ''.join(seqs)  # yield a fastq record
                    break

            if last:  # reach EOF before reading enough quality
                yield name, sequence, None  # yield a fasta record instead
                break


def write_fasta(counter, prefix, out_name):
    """
    Given a counter (typically containing collapsed reads and their abundances) write them out to a FASTA file.
    Reads will be output in order of abundance, with various fields to their identifier line:
     - A prefix (e.g. the donor file or conditions from which the sequences arose)
     - Rank (determined by order in the sorted counter
     - Abundance (determined by it's count post-collapsing)
    """

    with opener(out_name, 'w') as out_file:
        rank = 1
        for tr in counter.most_common():
            out_string = ">" + prefix + "_Rank:" + str(rank) + "_Abundance:" + str(tr[1]) \
                         + "_Frequency:" + str(tr[1] / sum(counter.values())) + "\n" + tr[0] + "\n"
            out_file.write(out_string)
            rank += 1


def get_qual_scores(qualstring):
    """ Convert FASTQ quality scores to their integer Q score equivalent """
    quality_list = [ord(x) - 33 for x in qualstring]
    return quality_list


def ldist(one, two):
    """ Return the Levenshtein distance of two strings """
    return lev.distance (one, two)


def collapse(fl):
    """
    Given a properly-formatted FASTQ file, return:
     - a counter containing the DNA sequenced stringently collapsed with their abundances, and
     - a counter containing information relating to the collapsing/error-correction
    Sequences are collapsed along their entire sequence, to strict criteria, including:
     - Each base in each barcode must be >= Q25
     - Each barcode must have at least three reads
     - The consensus base at each position is then determined
    NB: The first twelve nucleotides of each read much consist of random sequences incorporated before PCR
    """
    with opener(fl, 'r') as infile:
        counts = coll.Counter()
        barcodes = coll.defaultdict(list)
        ids = coll.defaultdict(list)

        for readid, seq, qual in readfq(infile):
            counts['total_seq_in'] += 1

            # Check no ambiguous base calls
            if "N" in seq:
                counts['N'] += 1
                continue

            # Read in barcode
            bc = seq[:12]

            # Check barcode is of sufficient quality
            bc_q_thresh = 25  # specify a phred score that each base call in the barcode must beat

            # Infer the location of the transcript information

            if all(b > bc_q_thresh for b in get_qual_scores(qual)[:12]):
                pad = 24
                transcript = seq[pad:-pad]

            else:
                counts['bc_fail'] += 1
                continue

            # Put successfully checked reads into dict based on their barcodes
            barcodes[bc].append(transcript)
            ids[bc].append(readid)
            counts['seqs_passed_to_barcodes'] += 1

    # Generate consensus sequences within barcode families, then collapse
    read_threshold = 3  # number of reads a barcode has to have to be collapsed

    cons_seqs = coll.Counter()
    cons_ids = coll.defaultdict(list)
    cons_bc = coll.defaultdict(list)

    for bc in barcodes:

        # Check that this barcode has the threshold number of reads
        if len(barcodes[bc]) < read_threshold:
            counts['bc_too_few_reads'] += 1
            continue

        # Read all sequences into a Counter & find the biggest (and if reads are all the same we can skip collapsing)
        seqs = coll.Counter()
        for x in barcodes[bc]:
            seqs[x] += 1

        if len(seqs) == 1:
            counts['bc_all_seq_same'] += 1
            cons_seqs[x] += len(barcodes[bc])
            cons_bc[x].append(bc)
            for t in ids[bc]:
                cons_ids[x].append(t)
            continue

        # If the most common sequence is a certain amount bigger than the next biggest clone, take that as genuine
        amount_bigger = 3
        top2 = seqs.most_common(2)

        if top2[0][1] > (top2[1][1] * amount_bigger):
            counts['bc_from_proto'] += 1
            cons_seqs[top2[0][0]] += len(barcodes[bc])
            cons_bc[top2[0][0]].append(bc)
            for t in ids[bc]:
                cons_ids[top2[0][0]].append(t)
            continue

        # Otherwise we'll have to generate a consensus the old fashioned way:
        # generate an array of reads and take the most common base at each position

        # First see whether all reads are the same length
        lens = coll.Counter()
        for s in seqs:
            lens[len(s)] += seqs[s]

        # Make sure that there are enough of same length to pass filter
        if lens.most_common(1)[0][1]:
            counts['bc_too_few_same_len'] += 1
            continue

        # Get consensus seq from any left
        arr_len = lens.most_common(1)[0][0]  # set the length of array as most frequently used length
        arr_height = lens.most_common(1)[0][1]  # set the length of array as most frequently used length

        # Make array: rows = # reads of same length, columns = # bases in reads
        read_arr = np.chararray((arr_height, arr_len))

        # Enter sequence characters into rows
        row = 0
        for sq in barcodes[bc]:
            if len(sq) == arr_len:
                read_arr[row:row + 1] = sq
                row += 1

        # Take consensus character and add to overall string
        consensus = ""
        for c in range(arr_len):
            base_count = coll.Counter(read_arr[:, c])
            # Check that most frequent base at a position is uniquely common (i.e. no draw)
            if len(base_count) > 1:
                if base_count.most_common(2)[0][1] == base_count.most_common(2)[1][1]:
                    base_count['N'] += base_count.most_common(1)[0][1] + 1
            # Otherwise, having passed the tests, add the most common base at that position to the consensus!
            consensus = consensus + base_count.most_common(1)[0][0]

        if 'N' in consensus:
            counts['bc_with_N'] += 1

        if len(consensus) == arr_len:
            cons_seqs[consensus] += len(barcodes[bc])
            cons_bc[consensus].append(bc)
            for t in ids[bc]:
                cons_ids[consensus].append(t)
            counts['bc_passed'] += 1
        else:
            counts['bc_wrong_len_cons'] += 1

    # Then collapse on barcodes
    dist_thresh = 3
    col_reads = coll.Counter()

    for cons in cons_bc.keys():
        tracking = coll.Counter()
        clustered = coll.Counter()
        cbc = cons_bc[cons]  # find all barcodes associated with a given consensus
        for b1 in cbc:
            if tracking[b1] == 0:
                ref_bc = b1
                tracking[b1] = 1
                clustered[b1] = 1
            for b2 in cbc:
                if tracking[b2] == 0:
                    if ldist(ref_bc, b2) < dist_thresh:
                        tracking[b2] = 1
                        clustered[ref_bc] += 1
        col_reads[cons] = len(clustered)

    return col_reads, counts


def dna_check(possible_dna):
    """
    :param possible_dna: A sequence that may or may not be a plausible DNA (translatable!) sequence
    :return: True/False
    """

    if possible_dna:
        return set(possible_dna.upper()).issubset({'A', 'C', 'G', 'T', 'N'})
    else:
        return False


def today():
    """
    :return: Today's date, in ISO format
    """
    return datetime.datetime.today().date().isoformat()


def out_folder(dir_suffix):
    """
    :return: The path to the plots directory subfolder for results plotted on this day (creating it if needed)
    """
    out_dir = output_dir + today() + '-' + dir_suffix + '/'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    return out_dir


def sub_out_folder(analysis_out_dir, dir_suffix):
    """
    :return: The path to the plots directory subfolder for results plotted on this day (creating it if needed)
    """
    if not analysis_out_dir.endswith('/'):
        analysis_out_dir += '/'
    out_dir = analysis_out_dir + dir_suffix + '/'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    return out_dir


def open_previous_data(directory_suffix, data_name, data_type):
    """
    :param directory_suffix: suffix of directory to look for the file in
    :param data_name: the named portion of a file produced by one of the upstream scripts
    :param data_type: file extension, e.g. 'tsv'
    :return: the most recent TSV file fitting that name (pd.df), plus the path to the root output directory
    """
    data_type = data_type.lower()
    if data_type.lower() not in ['tsv']:
        raise IOError('Inappropriate data type not in acceptable list (csv/pkl): \"' + data_type + '\"')

    all_dir_hits = [x for x in os.listdir(output_dir) if x.endswith(directory_suffix)]
    if all_dir_hits:
        all_dir_hits.sort()
        most_recent = output_dir + all_dir_hits[-1] + '/'
        print("\tFound the directory: " + most_recent)
        all_file_hits = [x for x in os.listdir(most_recent) if data_name.lower() in x.lower() and
                         (x.lower().endswith(data_type) or x.lower().endswith(data_type + '.gz'))]

        if len(all_file_hits) == 1:
            print("\t\tFound a single matching file: " + all_file_hits[0])
            data_file = all_file_hits[0]
        elif len(all_file_hits) == 0:
            raise IOError("Couldn't find any matching files for " + data_name + " in " + most_recent)
        elif len(all_file_hits) > 1:
            raise IOError("Found >1 matching files for " + data_name + " in " + most_recent)

        if data_type == 'tsv':
            return pd.read_csv(most_recent + data_file, compression = 'infer', sep='\t'), most_recent

        else:
            raise ValueError("Unknown data type: " + data_type)
    else:
        raise IOError('Cannot find a directory matching \"' + directory_suffix + '\" in ' + output_dir)


def run_bash(cmd):
    """
    Run a command via the bash terminal
    :param cmd: The command to be run
    :return: nothing
    """
    subprocess.call(cmd, shell=True)


codons = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
    'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
}

stitchr_headers = ['TCR_name', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3',
                   'TRAC', 'TRBC', 'TRA_leader', 'TRB_leader', 'Linker', 'Link_order',
                   'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq', 'TRB_3_prime_seq']

data_dir = '../Data/'
output_dir = '../Outputs/'

emerson_dir = data_dir + 'Emerson/'
raw_emerson_dir = emerson_dir + 'Raw/'
int_emerson_dir = emerson_dir + 'Intermediate/'
proc_emerson_dir = emerson_dir + 'Processed/'
bench_emerson_dir = emerson_dir + 'Benchmarking/'

vdjdb_dir = data_dir + 'VDJdb/'
raw_vdjdb_dir = vdjdb_dir + 'Raw/'
proc_vdjdb_dir = vdjdb_dir + 'Processed/'
bench_vdjdb_dir = vdjdb_dir + 'Benchmarking/'

heather_dir = data_dir + 'Heather/'
raw_heather_dir = heather_dir + 'Raw/'
int_heather_dir = heather_dir + 'Intermediate/'
proc_heather_dir = heather_dir + 'Processed/'
rerun_heather_dir = heather_dir + 'Processed-Rerun/'

immunesim_dir = data_dir + 'immuneSIM/'
raw_immunesim_dir = immunesim_dir + 'Raw/'
int_immunesim_dir = immunesim_dir + 'Intermediate/'
proc_immunesim_dir = immunesim_dir + 'Processed/'

supp_script_dir = 'Supplementary-Scripts/'

exts = ['png']  # If desired additional plot types can be produced by changing this value
sizes = [3.5, 4, 4.5, 5]  # Sizes of plots produced, for convenience when assembling figures

gits = {'stitchr': 'https://github.com/JamieHeather/stitchr.git',
        'immunoseq2airr': 'https://github.com/JamieHeather/immunoseq2airr.git',
        'autoDCR': 'https://github.com/JamieHeather/autoDCR.git'}

__all__ = [x for x in dir() if x.endswith('_dir')]
