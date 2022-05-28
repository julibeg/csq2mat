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


# %% ###################################################################
# parse the GFF
gff_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    "/test-data/refgenome/MTB-h37rv_asm19595v2-eg18.gff"
)
vcf_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    # "/test-data/amikacin_split_no_PE_PPE_filtered_csq.bcf.gz"
    "/test-data/test_subset_csq.vcf"
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
other_variants = []
missense_variants = {}


def index2bitmask(n: int) -> int:
    """
    Returns the integer corresponding to the bitmask used for the n-th consequence in
    the `INFO/BCSQ` field for sample having the homozygous genotype (i.e. `3` for the
    first consequence, `12` for the second, `48` for the third, etc.).
    """
    return 3 << (2 * n)


for _, entry in res.iterrows():
    if entry["CSQ"] == ".":
        var_id = "non_coding_" + "_".join(
            entry[["CHROM", "POS", "REF", "ALT"]].astype(str)
        )
        other_variants.append(
            pd.Series(
                entry[samples].apply(lambda x: x[0]).replace(".", 3),
                index=samples,
                name=var_id,
                dtype=pd.SparseDtype(np.int8, 0),
            )
        )
    else:
        # there might be more than one consequence in the CSQ field --> split them
        csqs = entry["CSQ"].split(",")
        for i, csq in enumerate(csqs):
            # skip if part of a compound consequence
            if csq.startswith("@"):
                continue
            # get the bitmask corresponding to this consequence (make sure that it's a
            # string; otherwise the comparison below will fail)
            bitmask = str(index2bitmask(i))
            csq, gene_name, *_, aa_change, _ = csq.split("|")
            if csq == "synonymous":
                var_id = f"synonymous_{gene_name}_{aa_change}"
                other_variants.append(
                    pd.Series(
                        entry[samples].apply(lambda x: x[0]).replace(".", 3),
                        index=samples,
                        name=var_id,
                        dtype=pd.SparseDtype(np.int8, 0),
                    )
                )
            elif csq == "missense":
                # create the variant ID
                ref_aa = aa_change.split(">")[0]
                var_id = f"missense_{gene_name}_{ref_aa}"
                # get the genotypes ('1' if the sample has the GT and the consequence
                # with the correct bitmask)
                 gts = entry[samples].apply(
                            lambda x: 3
                            if x[0] == "."
                            else 1
                            if x[0] == "1" and x[-1] == bitmask
                            else 0
                        )
                # check if we have already encountered a missense variant in the same
                # codon
                if var_id in missense_variants:
                    missense_variants[var_id] |= gts
                else:
                    missense_variants[var_id] = gts
                # THIS NEEDS TO BE TESTED PROPERLY!
            else:
                var_id = "other_" + "_".join(
                    entry[["CHROM", "POS", "REF", "ALT"]].astype(str)
                ) + '_' + csq
                other_variants.append(
                    pd.Series(
                        entry[samples].apply(lambda x: x[0]).replace(".", 3),
                        index=samples,
                        name=var_id,
                        dtype=pd.SparseDtype(np.int8, 0),
                    )
                )


final = pd.DataFrame(final)
final
# %% ###################################################################
bla = pd.Series([0, 1, 0])
blu = pd.Series([0, 1, pd.NA], dtype=pd.SparseDtype(float, 0))
bla | blu
# %% ###################################################################


def process_entry(row):
    if row["CSQ"] == ".":
        return row[samples].apply(lambda x: x[0])
    else:
        return row[samples].apply(lambda x: x[-1] if x[0] != "." else ".")


res2 = res.copy()
res2[samples] = res2.apply(process_entry, axis=1)
print(res[samples])
res2[samples]
