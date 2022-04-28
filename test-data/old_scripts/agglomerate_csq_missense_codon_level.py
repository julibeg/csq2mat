#!/usr/bin/env python3
"""
Takes a VCF or BCF file (must have been annotated with `bcftools csq`) and the
number of CPUs to use. Writes a dataframe (as CSV) to STDOUT with missense
variant IDs in the columns and sample IDs in the rows. When there was a missense
mutation in that codon, then there will be a '1' in the corresponding cell in
the dataframe. Otherwise there will be a '0'.
"""

import pandas as pd
import subprocess
import multiprocessing as mp
from tqdm import tqdm
import sys


def eprint(string):
    sys.stderr.write(string)
    sys.stderr.flush()


vcf_fname = sys.argv[1]
cpus = int(sys.argv[2]) if len(sys.argv) == 3 else 1

# get sample IDs
eprint('get sample IDs from variant file... ')
samples = pd.Series(subprocess.run(
    f'bcftools query -l {vcf_fname}'.split(), capture_output=True,
    universal_newlines=True).stdout.strip().split('\n'), name='samples')
eprint('done\n')

missense_df = pd.DataFrame(index=samples)


def query_and_process_VCF_chunk(vcf_fname, gene, length=0.9):
    """
    Function used for asynchronously querying the annotated variant file and
    returning codons with missense SNPs.
    """
    lof_samples = set()
    missense_dict = {}
    start, end = genes.loc[gene]
    # don't consider the last 10% of the gene
    gene_length = end - start
    new_end = start + int(round(0.9 * gene_length))
    query = subprocess.run(
        f"bcftools query -f '[%SAMPLE:%POS:%TBCSQ{{0}}\t]' {vcf_fname} "
        f"-r Chromosome:{start}-{end}".split(' '),
        capture_output=True, universal_newlines=True).stdout.strip().strip('\'')
    # eprint(f"bcftools query -f '[%SAMPLE:%POS:%TBCSQ{{0}}\t]' {vcf_fname} "
    #        f"-r Chromosome:{start}-{end}\n")
    # eprint(f'{len(query)}\n')
    for e in query.split('\t'):
        e = e.strip('\'')
        if e == '':
            continue
        s, pos, csqs = e.split(':')
        pos = int(pos)
        # in case of overlapping genes the consequences are separated by commas
        for csq in csqs.split(','):
            if csq.startswith('@'):
                continue
            if 'missense' in csq and not csq.startswith('*'):
                if csq.split('|')[0] != 'missense':
                    eprint(csq + '\n')
            csq, csq_gene = csq.split('|')[:2]
            if csq.startswith('*') or gene != csq_gene:
                # ignore mutations after a stop codon and we should only
                # consider the consequences of gene of interest
                continue
            if ('frameshift' in csq or 'stop_gained' in csq) and pos <= new_end:
                # handle LOF case
                lof_samples.add(s)
    return gene, lof_samples


# setup process pool and get LOF samples for every gene. then, set the
# samples to '1' in the df via `add_lof_samples` as callback. when querying
# the variant file, ignore the last 10% of every gene as LOF mutations in that
# region might not be as impactful.
eprint(f'query and process variants for each of {len(genes)} genes...\n')
with mp.Pool(cpus) as p, tqdm(total=len(genes)) as pbar:

    # define callback here to be able to use pbar
    def callback(output):
        gene, lof_samples = output
        lof_df.loc[lof_samples, gene] = 1
        pbar.update(1)

    for gene in genes.index:
        p.apply_async(func=query_and_process_gene,
                      args=(vcf_fname, gene),
                      callback=callback)
    p.close()
    p.join()

eprint('write results... ')
lof_df.to_csv(sys.stdout)
eprint('done\n')
