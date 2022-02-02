#!/usr/bin/env python3

import pysam
from collections import Counter
from pathlib import Path
import click

@click.command()
@click.argument('bam-file')
@click.option('--read-length', default=40, help='Minimum length of mapped reads to count in')
def count_full_mapped_reads(bam_file, read_length):
    bam = pysam.AlignmentFile(bam_file, "rb")
    full_reads = []
    for read in bam:
        if len(read.cigar) == 1:
            if (read.cigar[0][0] == 0 and read.cigar[0][1] >= read_length):
                full_reads.append(read.cigar[0][1])
    count_ = dict(Counter(full_reads))
    sum_mapped_reads = sum(count_.values())
    filename = Path(bam_file).name
    print(f"{filename}\t{sum_mapped_reads}")

if __name__ == "__main__":
    count_full_mapped_reads()
