import pandas as pd
import numpy as np
from typing import List, Tuple
import subprocess
import pathlib
import re
import io
import multiprocessing


def parse_gff(gff_fname: str) -> pd.DataFrame:
    """
    Parse a GFF file to get the start and end coordinates of all CDSs (assuming that
    there is only one CDS per gene entry)
    """
    genes = pd.DataFrame(columns=["name", "start", "end"])
    with open(gff_fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            # we are only interested in 'gene' and 'CDS' entries
            if (gff_elements := line.strip().split("\t"))[2] not in ("gene", "CDS"):
                continue
            _, _, typ, start, end, *_, desc = gff_elements
            desc_dict = dict(x.split("=") for x in desc.split(";"))
            if typ == "gene":
                # if the line represents a gene entry, add the gene ID and name to the
                # DataFrame
                gene_id = desc_dict["gene_id"]
                gene_name = desc_dict.get("Name", gene_id)
                genes.loc[gene_id, "name"] = gene_name
            elif typ == "CDS":
                # This is the CDS of the current gene entry --> add start and end to the
                # DataFrame
                prot_id = desc_dict["protein_id"]
                genes.loc[prot_id, ["start", "end"]] = int(start), int(end)
    genes[['start', 'end']] = genes[['start', 'end']].astype(int)
    return genes.sort_values("start")


def parse_VCF_header(vcf_file: str) -> Tuple[str, int, List[str]]
    """
    Parse the header of a bacterial VCF file to get the name and length of the
    chromosome as well as the sample IDs.
    """
    header = subprocess.run(
        ["bcftools", "view", "-h", vcf_file], capture_output=True, text=True
    ).stdout
    for line in header.strip().split("\n"):
        if "contig=" in line:
            # unpacking into a tuple of length 1 makes sure that there was only one
            # match
            (chr_name,) = re.findall("ID=(.*?),", line)
            (chr_length,) = re.findall("length=([0-9]+)", line)
            chr_length = int(chr_length)
    # the last line is the line holding the sample IDs
    samples = line.split("\t")[9:]
    return chr_name, chr_length, samples


