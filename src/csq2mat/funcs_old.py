import argparse
import pandas as pd
from typing import Tuple, List


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


def get_bitmask_value(n: int) -> int:
    """
    Returns the integer corresponding to the bitmask used for the n-th consequence in
    the `INFO/BCSQ` field for sample having the homozygous genotype (i.e. `3` for the
    first consequence, `12` for the second, `48` for the third, etc.).
    """
    return 3 << (2 * n)


def parse_GFF(gff_file: str) -> None:
    # -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a gff file and return two DataFrames with genomic start and end positions;
    one for the genes found in the GFF and one for the intergenic regions.

    Parameters
    ----------
    gff_file : str
        A GFF file.

    Returns
    -------
    gene_positions : pd.DataFrame
        Start and end positions of the genes found in the GFF.
    intergenic_positions : pd.DataFrame
        Start and end positions of the intergenic regions inferred from the GFF.
    """
    with open(gff_file, "r") as gff:
        # parse the header first
        header: List[str] = []
        for line in gff:
            if line.startswith("#"):
                header.append(line)
            else:
                break
        # make sure the file is GFF version 3
        if "gff-version 3" not in header[0].lower():
            raise ValueError(
                "The GFF appears not to be GFF version 3. If it is, please add "
                'the line "##gff-version 3" to the top of the file.'
            )
        # now parse the remaining lines and create the lists for the genes and
        # intergenic regions
        prev_gene: str = ''
        for line in gff:
            if line.startswith('#'):
                continue
            _, _, entry_type, start, end, *desc = line.split()
            if entry_type is gene:
                
                