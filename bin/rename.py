#!/usr/bin/env python3

import pandas as pd
# import gspread
import sys
import os
import subprocess
from pathlib import Path
import re

def rename(lookup_file, fastq):
    """
    Rename barcode by cog number
    """
    lktbl = None
    if lookup_file.endswith(".csv"):
        lktbl = pd.read_csv(lookup_file)
    elif lookup_file.endswith(".xlsx"):
        lktbl = pd.read_excel(lookup_file)
    else:
        sys.exit("Lookup table file has incorrect format.\ Should be csv/tsv or excel!")
    
    fastq_name = re.findall(r'barcode[0-9]{2,4}', fastq)
    if len(fastq_name) > 0:
        lktbl = lktbl.applymap(lambda x: x.strip().replace(" ", "_") if isinstance(x, str) else x)
        lktbl["barcode"] = lktbl["barcode"].str.lower()
        lktbl["barcode"] = lktbl["barcode"].str.replace("bc", "barcode")

        cog_value = lktbl[lktbl["barcode"] == fastq_name[0]]["cog"].to_list()
        if len(cog_value) == 1:
            cwd = Path.cwd()
            src = cwd.joinpath(fastq)
            src.rename("nanopore.fq")
            dst_name = f"{cog_value[0]}.fastq"
            Path(dst_name).symlink_to("nanopore.fq")
            # src.unlink()
    else:
        os.exit(0)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        lookup_file = sys.argv[1]
        fastq = sys.argv[2]
        rename(lookup_file, fastq)
    else:
        sys.exit("Wrong input!")
