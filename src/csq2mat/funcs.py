# %% ###################################################################
import argparse
import pandas as pd
import numpy as np
import subprocess
from typing import Tuple, List
import functools


def get_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        csq2mat: A tool for LOF burden testing and combining missense mutations on the
        codon level for bacterial GWAS and phenotype prediction.
        """,
    )
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        required=True,
        help="VCF/BCF file (must have been annotated with `bcftools csq`) [required]",
        metavar="FILE",
    ),
    parser.add_argument(
        "-g",
        "--gff",
        type=str,
        required=True,
        help=(
            "GFF file (must be the same file that was used to annotate "
            "the VCF) [required]"
        ),
        metavar="FILE",
    ),
    return parser.parse_args()


def index2bitmask(n: int) -> int:
    """
    Returns the integer corresponding to the bitmask used for the n-th consequence in
    the `INFO/BCSQ` field for sample having the homozygous genotype (i.e. `3` for the
    first consequence, `12` for the second, `48` for the third, etc.).
    """
    return 3 << (2 * n)


def get_leftmost_set_bit(n: int) -> int:
    if n < 0:
        raise (
            ValueError,
            f"Cannot get leftmost set bit for negative number: number was {n}",
        )
    elif n == 0:
        return 0
    lmb = 1
    while n > 1:
        n >>= 1
        lmb += 1
    return lmb


@functools.lru_cache(maxsize=None)
def bitmask2index(bm: int) -> int:
    # make sure the bitmask was `3` or greater
    if bm < 3:
        raise ValueError(
            f"Cannot find an index for a bitmask smaller than 3: bitmask was {bm}"
        )
    # get the leftmost set bit first
    lmb = get_leftmost_set_bit(bm)
    # make sure the bitmask was valid (i.e. two consecutive bits set due to homozygosity
    # in the VCF and the position of the leftmost set bit was a multiple of 2)
    reconstructed_bm = (1 << (lmb - 1)) + (1 << (lmb - 2))
    if (reconstructed_bm != bm) or lmb % 2:
        raise ValueError(f"Looks like an invalid bitmask: bitmask was {bm}")
    return lmb // 2


def get_sample_ids(VCF_file: str) -> List[str]:
    p = subprocess.run(
        ["bcftools", "query", VCF_file, "-l"],
        capture_output=True,
        text=True,
    )
    if p.returncode:
        raise RuntimeError(f'ERROR getting sample IDs from "{VCF_file}":\n{p.stderr}')
    return p.stdout.strip().split()


def query_VCF(VCF_file: str, query: str, *extra_args) -> str:
    p = subprocess.run(
        ["bcftools", "query", VCF_file, "-f", query, *extra_args],
        capture_output=True,
        text=True,
    )
    if p.returncode:
        raise RuntimeError(
            f'ERROR querying VCF file "{VCF_file}" with query "{query}":\n{p.stderr}'
        )
    return p.stdout


def parse_genotypes(gt_str: str):
    for gt in gt_str.strip(",").split(","):
        if (g := gt[0]) == ".":
            yield np.nan
        else:
            yield int(g)


def parse_VCF_query(query: str):
    for line in query.strip().split("\n"):
        pos, ref, alt, csqs, gts, csq_bitmasks = line.strip().split("\t")
        varID = "_".join((pos, ref, alt))
        if csqs == ".":
            # variant is not in a CDS --> simply return the parsed genotypes
            yield pd.Series(
                parse_genotypes(gts), name=varID, dtype=pd.SparseDtype(float, 0)
            )
        # csqs = csqs.split(",")
        # for csq in csqs:
        #     if


test_vcf = (
    "/home/julsy/git/csq2mat/test-data/"
    "amikacin_split_no_PE_PPE_filtered_csq_top10kb.bcf.gz"
)
test_query = "%POS\t%REF\t%ALT\t%INFO/BCSQ\t[%GT,]\t[%BCSQ,]\n"
q = query_VCF(test_vcf, test_query)
res_df = pd.concat(parse_VCF_query(q), axis=1)
res_df.index = get_sample_ids(test_vcf)
res_df
# %% ###################################################################
import numpy as np

df = pd.DataFrame(np.zeros((4, 3)), dtype=pd.SparseDtype(float, 0))
df[3] = pd.Series([1, np.nan, 0, 0], dtype=pd.SparseDtype(float, 0))
df = pd.concat((df, pd.DataFrame([0, 1, 0, 0]).T), axis=0)
df.dtypes
# %% ###################################################################
df = pd.DataFrame(np.zeros((4, 3)))
df = pd.concat((df, pd.Series([1, 0, 0])), axis=0)
df
# %% ###################################################################
pd.Series(parse_genotypes("0/0,1/1,./."), dtype=pd.SparseDtype(float, 0), name="bla")
# %% ###################################################################
