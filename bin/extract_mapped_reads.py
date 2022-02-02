#!/usr/bin/env python3
import pysam
import re
import logging
import click

# logging.setLevel(logging.INFO)


def extract_mapped_reads(bamfile):
    infile = pysam.AlignmentFile(bamfile, "rb")
    for segment in infile.fetch(until_eof=True):
        print(
            f"{segment.query_name}\t{segment.query_alignment_length}\t1")


@click.command()
@click.argument('bamfile')
def main(bamfile):
    extract_mapped_reads(bamfile)


if __name__ == "__main__":
    main()
