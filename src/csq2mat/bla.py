# type: ignore
# %% ###################################################################
import pandas as pd
import numpy as np
import subprocess
import pathlib
import re
import io
import multiprocessing
from tqdm import tqdm
import funcs


LOF_FRAMESHIFT_STOP_GAINED_THRESHOLD = 0.9
LOF_INFRAME_DELETION_LENGTH_TRHESHOLD = 30
LOF_INFRAME_INSERTION_LENGTH_TRHESHOLD = 30


# %% ###################################################################
# parse the GFF
gff_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    "/test-data/refgenome/MTB-h37rv_asm19595v2-eg18.gff"
)
vcf_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    # "/test-data/amikacin_split_no_PE_PPE_filtered_csq.bcf.gz"
    # "/test-data/test_subset_csq.vcf"
    "/test-data/toy_csq.vcf"
)
# %% ###################################################################
# parse the VCF header
chr_name, chr_length, samples = funcs.parse_VCF_header(vcf_file)
# %% ###################################################################
query = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n"
res = pd.read_csv(
    io.StringIO(funcs.run_bcftools_cmd(["query", "-f", query], vcf_file)),
    sep="\t",
    header=None,
    names=["CHROM", "POS", "REF", "ALT", "CSQ", *samples],
)
print(res)
# %% ###################################################################


def csq_index_to_bitmask(n: int) -> int:
    """
    Returns the integer corresponding to the bitmask used for the n-th consequence in
    the `INFO/BCSQ` field for sample having the homozygous genotype (i.e. `3` for the
    first consequence, `12` for the second, `48` for the third, etc.).
    """
    return 3 << (2 * n)


def csq_bitmasks() -> int:
    """
    A generator yielding INFO/BCSQ bitmasks in increasing order for samples with
    homozygous genotypes (i.e. `3`, `12`, `48`, etc.).
    """
    n = 0
    while True:
        yield 3 << (2 * n)
        n += 1


def is_LOF(csq: str, aa_change: str) -> bool:
    """
    Check if the variant is a loss of function variant according to our definitions.
    """
    # if the variant is after the stop codon --> no LOF
    if csq.startswith("*"):
        return False
    # frameshift and early stop codons
    if "frameshift" in csq or "stop_gained" in csq:
        # determine if the mutation was early enough in the gene to be considered LOF
        codon_num, orig_gene_aa_length = map(
            int, re.match(r"(\d+)[a-zA-Z]\.\.(\d+)>", aa_change).group(1, 2)
        )
        if codon_num < LOF_FRAMESHIFT_STOP_GAINED_THRESHOLD * orig_gene_aa_length:
            return True
    # large deletions
    if "inframe_deletion" in csq:
        # if there is no '..' in aa_change, then it was a small deletion --> not LOF
        if ".." in aa_change:
            # determine the length of the deletion in number of amino acids
            del_start, del_end = map(
                int, re.match(r"(\d+)[a-zA-Z]\.\.(\d+)>", aa_change).group(1, 2)
            )
            deletion_length = del_end - del_start - 1
            if deletion_length > LOF_INFRAME_DELETION_LENGTH_TRHESHOLD:
                return True
    # large insertions
    if "inframe_insertion" in csq:
        # if there is no '..' in aa_change, then it was a small insertion --> not LOF
        if ".." in aa_change:
            # determine the length of the insertion in number of amino acids
            ins_start, ins_end = map(
                int, re.search(r">(\d+)[a-zA-Z]\.\.(\d+)", aa_change).group(1, 2)
            )
            insertion_length = ins_end - ins_start - 1
            if insertion_length > LOF_INFRAME_INSERTION_LENGTH_TRHESHOLD:
                return True
    return False


# %% ###################################################################
other_variants = []
missense_variants = {}
LOF_variants = {}

for _, entry in res.iterrows():
    # there might be more than one consequence in the CSQ field --> split them and
    # handle each consequence individually
    csqs = entry["CSQ"].split(",")
    # loop over the consequences and the corresponding bitmasks (i.e. bitmasks for
    # samples with that consequence)
    for csq, bitmask in zip(csqs, csq_bitmasks()):
        # skip if part of a compound consequence (the compound consequence will be
        # handled at some point anyway)
        if csq.startswith("@"):
            continue
        # get a vector with `1`s for the samples having this particular consequence
        has_csq = (
            entry[samples]
            .apply(
                # we use `3` to represent non-calls --> this also works with the binary
                # `OR` operation in the missense section below (`3 | 0` and `3 | 1` will
                # both give `3`), but must be replaced by NaN at the end.
                lambda x: 3
                if x[0] == "."
                else 1
                if int(x.split(":")[1]) == bitmask
                else 0
            )
            .astype(pd.SparseDtype(int, 0))
        )
        # split the CSQ field
        csq, gene_name, *_, aa_change, _ = csq.split("|")
        # now handle the different types of consequences. For non-coding, synonymous,
        # and non-missense CDS variants the process will be the same and only the
        # variant ID in the output will change. We'll start with the LOF consequences.
        if is_LOF(csq, aa_change):
            var_id = f"LOF_{gene_name}"
            if var_id in LOF_variants:
                LOF_variants[var_id] |= has_csq
            else:
                LOF_variants[var_id] = has_csq
        elif csq == "missense":
            # now process the missense consequences
            ref_aa = aa_change.split(">")[0]
            var_id = f"missense_{gene_name}_{ref_aa}"
            # check if we have already encountered a missense variant in the same
            # codon --> if so, combine with the current variant
            if var_id in missense_variants:
                missense_variants[var_id] |= has_csq
            else:
                missense_variants[var_id] = has_csq
        else:
            # now, let's handle the rest (synonymous, non-coding, and non-missense CDS
            # variants)
            if csq == "synonymous":
                var_id = f"synonymous_{gene_name}_{aa_change}"
            elif csq == '.':
                var_id = "non_coding_" + "_".join(
                    entry[["CHROM", "POS", "REF", "ALT"]].astype(str)
                )
            else:
                var_id = (
                    "other_"
                    + "_".join(entry[["CHROM", "POS", "REF", "ALT"]].astype(str))
                    + "_"
                    + csq
                )
            has_csq.name = var_id
            other_variants.append(has_csq)


print(res[["CSQ"] + samples])
final = pd.concat(
    [*other_variants, pd.DataFrame(LOF_variants), pd.DataFrame(missense_variants)],
    axis=1,
).T
final

