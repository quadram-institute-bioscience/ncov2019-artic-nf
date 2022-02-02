#!/usr/bin/env python
from collections import Counter
import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import subprocess
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import click
from loguru import logger


def count_full_mapped_reads(bam_file, read_length=148):
    bam = pysam.AlignmentFile(bam_file, "rb")
    sum_mapped_reads = bam.mapped
    full_reads = []
    long_read = False
    for read in bam:
        # logger.info(
        #     f"{read.infer_query_length()}-{read.infer_read_length()}-{read.query_alignment_length}-{read.cigar[0][1]}")
        if read.infer_query_length() > 250:
            long_read = True
            break
        if (read.query_alignment_length >= read_length):
            full_reads.append(read.query_alignment_length)
    
    if not long_read:
        count_ = dict(Counter(full_reads))
        sum_mapped_reads = sum(count_.values())
    
    logger.info(f"Mapped reads: {sum_mapped_reads}")
    return sum_mapped_reads

@click.command()
@click.argument('bam')
def go(bam):
    logger.info("Starting...")
    count_full_mapped_reads(bam)

if __name__ == "__main__":
    go()
