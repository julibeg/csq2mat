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

# %% ###################################################################
# parse the GFF
gff_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    "/test-data/refgenome/MTB-h37rv_asm19595v2-eg18.gff"
)
vcf_file = (
    f"{pathlib.Path(__file__).resolve().parent.parent.parent}"
    # "/test-data/amikacin_split_no_PE_PPE_filtered_csq.bcf.gz"
    "/test-data/test_subset_csq.bcf.gz"
)
# %% ###################################################################
# parse the GFF to get the start and end coordinates of all protein-coding genes
genes = pd.DataFrame(columns=["name", "start", "end"])
with open(gff_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        # we are only interested in 'gene' and 'CDS' entries
        if (line := line.strip().split("\t"))[2] not in ("gene", "CDS"):
            continue
        _, _, typ, start, end, *_, desc = line
        desc = dict(x.split("=") for x in desc.split(";"))
        if typ == "gene":
            # if the line represents a gene entry, add the gene ID and name to the
            # DataFrame
            gene_id = desc["gene_id"]
            gene_name = desc.get("Name", gene_id)
            genes.loc[gene_id, "name"] = gene_name
        elif typ == "CDS":
            # This is the CDS of the current gene entry --> add start and end to the
            # DataFrame
            prot_id = desc["protein_id"]
            genes.loc[prot_id, ["start", "end"]] = int(start), int(end)
genes.sort_values("start", inplace=True)
genes
# %% ###################################################################
# parse VCF header to get the name of the chromosome and the sample IDs
header = subprocess.run(
    ["bcftools", "view", "-h", vcf_file], capture_output=True, text=True
).stdout
for line in header.strip().split("\n"):
    if "contig=" in line:
        # unpacking into a tuple of length 1 makes sure that there was only one match
        (chr_name,) = re.findall("ID=(.*?),", line)
        (chr_length,) = re.findall("length=([0-9]+)", line)
        chr_length = int(chr_length)
# the last line is the line holding the sample IDs
samples = line.split("\t")[9:]
# %% ###################################################################
# find all intergenic regions (intervals not covered by any protein-coding gene)


def find_intergenic_regions(genes):
    old_end = 1
    old_gene_id = f"{chr_name}_Start"
    intergenic_regions = []
    for gene_id, gene in genes.iterrows():
        if gene["start"] > old_end:
            # the intervals taken by `bcftools -r` are inclusive --> adjust the
            # coordinates of the intergenic regions accordingly (i.e. +1 and -1 for
            # start and end)
            intergenic_regions.append(
                (
                    old_gene_id,
                    gene_id,
                    old_end + 1,
                    gene["start"] - 1,
                )
            )
        old_gene_id = gene_id
        old_end = max(old_end, gene["end"])
    if chr_length > old_end:
        intergenic_regions.append((gene_id, f"{chr_name}_End", old_end + 1, chr_length))
    return intergenic_regions


intergenic_regions = (
    pd.DataFrame(
        find_intergenic_regions(genes),
        columns=["upstream_gene", "downstream_gene", "start", "end"],
    )
    .set_index(["upstream_gene", "downstream_gene"])
    .astype(int)
)
intergenic_regions
# %% ###################################################################
# write genes + intergenic regions to CSV for checking them with bedtools
# genes.to_csv('check-intergenic-regions/genes.csv')
# intergenic_regions.to_csv('check-intergenic-regions/intergenic_regions.csv')
# %% ###################################################################
# process the intergenic regions


class BCFtoolsError(Exception):
    def __init__(self, msg):
        super().__init__(msg)

    @classmethod
    def from_CalledProcessError(cls, error):
        msg = f"{error.__str__()}\nSTDOUT: {error.stdout}\nSTDERR: {error.stderr}"
        return cls(msg)


def process_intergenic(start, end):
    query_fmt_str = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT]\n"
    region_str = f"{chr_name}:{start}-{end}"
    try:
        p = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                query_fmt_str,
                vcf_file,
                "-r",
                region_str,
                "--regions-overlap",
                "pos",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise BCFtoolsError.from_CalledProcessError(e)
    res = pd.read_csv(
        io.StringIO(p.stdout.strip()),
        sep="\t",
        header=None,
        names=["chr", "pos", "ref", "alt", "csq"] + samples,
        index_col=[0, 1, 2, 3],
    )
    if (
        not res["csq"]
        .apply(
            lambda x: x == "."
            or "non_coding" in x
            or "coding_sequence" in x
        )
        .all()
    ):
        # TODO: instead of permitting "non_coding" and "coding_sequence" consequences in
        # the intergenic regions, we could write all non-"." CSQs out to a logfile (and
        # not throw any error at all)
        bad_vars = res.query('csq != "."').index
        raise BCFtoolsError(
            f"INFO/CSQ in supposed intergenic region ({region_str}) was not empty! "
            + "The variants in question were: "
            + f'{list(("_".join(str(x) for x in idx) for idx in bad_vars))}.'
        )
    res = (
        res.drop(columns="csq")
        .applymap(lambda x: float(x[0]) if x[0] != "." else np.nan)
        .astype(pd.SparseDtype(float, 0))
    )
    return res


with multiprocessing.Pool(10) as pool, tqdm(total=len(intergenic_regions)) as pbar:
    results = [
        pool.apply_async(process_intergenic, args=(start, end), callback=pbar.update(1))
        for start, end in intergenic_regions.values
    ]
    results = pd.concat(r.get() for r in results)
results
# %% ###################################################################
bla = pd.concat(results)
# %% ###################################################################
pos = 1499211
intergenic_regions.query("start <= @pos and end >= @pos")
genes.loc[["Rv1329c", "Rv1330c"]]
# %% ###################################################################
pos = 1499211
genes.index.get_loc(genes.query("start <= @pos and end >= @pos").index[0])
genes.iloc[235:240]
# %% ###################################################################
find_intergenic_regions(genes.iloc[194:198])
# %% ###################################################################
start, end = intergenic_regions[0]
query_fmt_str = "%POS\t%REF\t%ALT\t%INFO/BCSQ\t[%GT,]\n"
region_str = f"{chr_name}:{start}-{end}"
try:
    p = subprocess.run(
        [
            "bcftools",
            "query",
            "-f",
            query_fmt_str,
            vcf_file,
            "-rx",
            region_str,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
except subprocess.CalledProcessError as e:
    bla = e
    print(dir(e))
    pass
# print(p.stdout.strip())
res = pd.read_csv(
    io.StringIO(p.stdout.strip()),
    sep="\t",
    names=["pos", "ref", "alt", "csq", "gt"],
    index_col=[0, 1, 2],
)
res["csq"] = "bla"
if not (res["csq"] == ".").all():
    bad_vars = res.query('csq != "."').index
    raise ValueError(
        f"Found string in INFO/CSQ in supposed intergenic region ({region_str}):"
        + "The variants in question were:\n"
        + "\n".join(",".join(str(x) for x in idx) for idx in bad_vars)
    )
