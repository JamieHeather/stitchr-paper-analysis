# -*- coding: utf-8 -*-

"""
alleles.py

Perform rudimentary novel allele inference on the Heather et al. TCRseq dataset

"""

import os
import collections as coll
import functions as fxn
import sys
sys.path.insert(0, fxn.supp_script_dir + 'autoDCR')
import inferTCR as infer


__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.0'
__author__ = 'Jamie Heather'

if __name__ == "__main__":
    scripts_dir = fxn.check_scripts_cwd()

    # First run autoDCR on each individual donor collapsed TCR repertoire
    dcr_dir = fxn.supp_script_dir + 'autoDCR/'
    rep_files = [x for x in os.listdir(fxn.int_heather_dir) if x.startswith('HV') and '.fasta' in x
                 and 'merged' not in x]
    rep_files.sort()

    for fasta in rep_files:
        # Run autoDCR with allele detection mode enabled ...
        fxn.run_bash('python3 ' + dcr_dir + 'autoDCR.py -fq ' + fxn.int_heather_dir + fasta \
            + ' -o ' + fxn.int_heather_dir + ' -dd ' + dcr_dir + ' -or forward -ad -jv')

        # ... and then infer potential novel TCRs
        ad_file = fasta.split('.')[0]
        fxn.run_bash('python3 ' + dcr_dir + 'inferTCR.py -in ' +
                     fxn.int_heather_dir + ad_file + '_infer-alleles.tsv.gz' + ' -dd ' + dcr_dir)

    # Compile the potential inferred alleles
    data_files = [x for x in os.listdir(os.getcwd()) if
                  x.startswith('HV') and x.endswith('infer-alleles_alleles.fasta')]
    data_files.sort()

    novel = coll.Counter()
    novel_lst = coll.defaultdict(list)
    seqs = coll.defaultdict()

    for df in data_files:
        with open(df, 'r') as in_file:
            for read_id, seq, qual in fxn.readfq(in_file):
                novel[read_id] += 1
                novel_lst[read_id].append(df)
                if read_id not in seqs:
                    seqs[read_id] = seq
                elif seqs[read_id] != seq:
                    raise ValueError("Mismatch between identified novel alleles with the same name!")

    with open(dcr_dir + 'compiled-inferred.fasta', 'w') as out_file:
        for gene in seqs:
            out_file.write(infer.imgt_fastafy(gene, seqs[gene], False))

    # Tidy up
    fxn.run_bash('tar cvfz allele-inference-files.tar.gz HV*; '
                 'rm HV*; '
                 'mv allele-inference-files.tar.gz ' + fxn.int_heather_dir)

    # Generate 'human-plus' autoDCR reference
    plus = 'human-plus.fasta'
    if 'compiled-inferred.fasta' in os.listdir(dcr_dir):
        fxn.run_bash('cat ' + dcr_dir + 'human.fasta ' + dcr_dir + 'compiled-inferred.fasta > ' +
                     dcr_dir + plus)
        fxn.run_bash('python3 ' + dcr_dir + 'generate-tag-files.py -in ' + dcr_dir + plus)

    else:
        raise IOError("'compiled-inferred.fasta' not in " + dcr_dir + "\nCode cannot proceed.")

    # Also add to the stitchr additional-genes file
    fxn.run_bash('cat ' + dcr_dir + 'compiled-inferred.fasta >> ' +
                 fxn.supp_script_dir + 'stitchr/Data/additional-genes.fasta')

    # Then the replotting.py script will use these data in the next step
