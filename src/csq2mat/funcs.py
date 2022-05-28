from __future__ import annotations  # for type hints of return value of constructor
import pandas as pd
import numpy as np
from typing import List, Tuple
import subprocess
import re


class BCFtoolsError(Exception):
    def __init__(self, msg: str) -> None:
        super().__init__(msg)

    @classmethod
    def from_CalledProcessError(
        cls, error: subprocess.CalledProcessError
    ) -> BCFtoolsError:
        msg = f"{error.__str__()}\nSTDOUT: {error.stdout}\nSTDERR: {error.stderr}"
        return cls(msg)


def run_bcftools_cmd(cmd: List[str], vcf_file: str) -> str:
    """
    Run a bcftools command on a VCF file.
    """
    try:
        return subprocess.run(
            ["bcftools", *cmd, vcf_file], capture_output=True, text=True
        ).stdout
    except subprocess.CalledProcessError as e:
        raise BCFtoolsError.from_CalledProcessError(e)
    

def parse_VCF_header(vcf_file: str) -> Tuple[str, int, List[str]]:
    """
    Parse the header of a bacterial VCF file to get the name and length of the
    chromosome as well as the sample IDs.
    """
    header = run_bcftools_cmd(["view", "-h"], vcf_file)
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


# def query_VCF(vcf_file):
#     """
#     Get the variant data, consequences, and genotypes from a VCF file.
#     """
#     query = "%CHROM\tPOS\tID\tREF\tALT\tINFO/BCSQ[\t%GT:BCSQ]\n"
#     res = run_bcftools_cmd(["query", "-f", query], vcf_file)
