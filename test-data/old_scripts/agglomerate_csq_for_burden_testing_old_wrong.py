#!/usr/bin/env python3
"""
Takes a VCF or BCF file (must have been annotated with `bcftools csq`), a
GFF file (ideally the one that has been used for the annotation), and the
number of CPUs to use. Writes a dataframe (as CSV) to STDOUT with genes
in the columns and samples in the rows. When a loss-of-function mutation
occurred in the first 90% of a gene in a particular sample, there will be a '1'
in the corresponding cell in the dataframe. Otherwise there will be a '0'.
"""

import pandas as pd
import subprocess
import multiprocessing as mp
from tqdm import tqdm
import sys


def eprint(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def get_LOF_samples(vcf_fname, start, end):
    """
    Function used for asynchronous querying the annotated variant file.
    """
    samples = set()
    query = subprocess.run(
        f"bcftools query -f '[%SAMPLE:%TBCSQ{{0}}\t]' {vcf_fname} "
        f"-r Chromosome:{start}-{end}".split(' '),
        capture_output=True, universal_newlines=True).stdout.strip()
    for e in query.split('\t'):
        e = e.strip('\'')
        if e == '':
            continue
        s, csq = e.split(':')
        csq = csq.split('|')[0]
        if csq.startswith('*'):
            # ignore variants downstream of stop codons
            continue
        if 'frameshift' in csq or 'stop_gained' in csq:
            samples.add(s)
    return samples


vcf_fname = sys.argv[1]
gff_fname = sys.argv[2]
cpus = int(sys.argv[3]) if len(sys.argv) == 4 else 1

# get sample IDs
eprint('get sample IDs from variant file... ')
samples = pd.Series(subprocess.run(
    f'bcftools query -l {vcf_fname}'.split(), capture_output=True,
    universal_newlines=True).stdout.strip().split('\n'), name='samples')
eprint('done\n')

# parse gff
eprint('parse GFF... ')
genes = pd.DataFrame(columns=['start', 'end'], dtype=int)
genes.index.name = 'genes'
with open(gff_fname, 'r') as gff:
    for line in gff:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        obj = line[2]
        if obj != 'gene':
            continue
        start, end = [int(x) for x in line[3:5]]
        notes = dict(map(lambda x: x.split('='), line[-1].split(';')))
        try:
            gene = notes['Name']
        except KeyError:
            gene = notes['ID'].split('gene:')[1]
        genes.loc[gene] = start, end
# setup result df
lof_df = pd.DataFrame(index=samples, columns=genes.index, data=0)
eprint('done\n')


# setup process pool and get LOF samples for every gene. then, set the
# samples to '1' in the df via `add_lof_samples` as callback. when querying
# the variant file, ignore the last 10% of every gene as LOF mutations in that
# region might not be as impactful.
eprint(f'query and process variants for each of {len(genes)} genes...\n')
with mp.Pool(cpus) as p, tqdm(total=len(genes)) as pbar:
    for gene, (start, end) in genes.iterrows():
        # define the callback function to add the results to the df within the
        # for loop and with the current gene as a default argument. this
        # makes sure that the correct value in the gene variable is used
        # by the callbacks. lambda expressions would use the final gene.
        def callback(lof_samples, gene=gene):
            lof_df.loc[lof_samples, gene] = 1
            pbar.update(1)
        # don't consider the last 10% of the gene
        gene_length = end - start
        new_end = start + int(round(0.9 * gene_length))
        p.apply_async(func=get_LOF_samples,
                      args=(vcf_fname, start, new_end),
                      callback=callback)
    p.close()
    p.join()

eprint('write results... ')
lof_df.to_csv(sys.stdout)
eprint('done\n')
