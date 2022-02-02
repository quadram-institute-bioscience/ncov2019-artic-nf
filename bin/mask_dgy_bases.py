#!/usr/bin/env python
# Thanh Le Viet - Quadram Institute Bioscience
# This script is used to trim off the dodgy bases which cause ambiguous nuc or artefact mutations in consensus sequence
# as they are not real, double checked in Positive controls between V3 and V4 runs
#
# 28/10/2021: Initial version


import click
import pysam
import pathlib
from collections import defaultdict
from loguru import logger
import os
from mycigar import Cigar

# Dodgy positions
POS = {
    # '8829': ["C", "A"],
    '8835': ["T", "C"],
    # '13791': ["A", "C"],
    '15521': ["T", "A"],
    # '21987': ["G", "A"],
    # '22227': ["C", "T"],
    # '23557': ["C", "T"],
    # '24410': ["G", "A"],
    # '25337': ["G", "C"],
    # '29870': ["C", "A"],
}


@click.command()
@click.argument('bam_file')
@click.option('-o', '--output-name', default=None, show_choices=True)
@click.option('--verbose', is_flag=True, default=False, show_default=True)
@click.option('-l', '--cutoff-length', default=20, help='Cut-off length of the edge containing the dodge base', show_default=True)
@click.option('--ref-name', default="MN908947.3", help='NCBI Reference name', show_default=True)
@click.option('-t', '--threads', default=4, help='# threads', show_default=True)
def soft_mask_dodgy_bases(bam_file, verbose, cutoff_length, ref_name, output_name, threads):
    if verbose:
        logger.enable("__main__")
    if output_name is None:
        output_name = pathlib.Path(bam_file).name
    if not output_name.endswith(".bam"):
        output_name = f"{output_name}.bam"
    assert pathlib.Path(bam_file).exists()

    logfile = f"{output_name}_trim_dodgy_bases.log"
    if pathlib.Path(logfile).exists():
        pathlib.Path(logfile).unlink()
    logger.add(logfile)

    os.system(f"samtools sort -@ {threads} -o trimmed.sorted.bam {bam_file}")
    os.system("samtools index trimmed.sorted.bam")

    logger.info(f"Bam file: {bam_file}")

    #1st query of reads in the regions of interest
    bam = pysam.AlignmentFile("trimmed.sorted.bam", "rb")
    # pysam.IndexedReads(bam)
    qseqs = defaultdict()
    for locus in POS:
        # logger.info(f"POS: {locus} - {POS[locus]} - {POS[locus][1]}")
        pos_of_interest = int(locus)
        # logger.info(pos_of_interest)
        for pileupcolumn in bam.pileup(ref_name, start=pos_of_interest-1, end=pos_of_interest+1, truncate=True, ignore_overlaps=False, stepper="nofilter"):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupcolumn.pos == pos_of_interest - 1:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == POS[locus][1]:
                            qseqs[pileupread.alignment.qname] = {'poi': pos_of_interest,
                                                                 'cigarstring': pileupread.alignment.cigarstring
                                                                 }
    bam.close()

    bam = pysam.AlignmentFile("trimmed.sorted.bam", "rb")
    fixed_bam = pysam.AlignmentFile("tmp.bam", "wb", template=bam)
    found_reads = list(qseqs.keys())
    n_fixed_reads = defaultdict()
    for read in bam:
        new_read = read
        if new_read.qname in found_reads and new_read.reference_start <= qseqs[new_read.qname]['poi'] - 1 <= new_read.reference_end:
            # logger.info(
            # f"{read.qname} mapped positions: [{read.reference_start} - {read.reference_end}] - Reverse: {read.is_reverse} - POI: {qseqs[read.qname]['poi']} - Soft-clip length: {qseqs[read.qname]['poi'] -  read.reference_start} CIGAR:{read.cigarstring} CIGAR tuples: {read.cigartuples}")
            soft_clip_length = qseqs[new_read.qname]['poi'] - \
                new_read.reference_start

            if soft_clip_length < cutoff_length:
                cigar_string = Cigar(new_read.cigarstring)
                cigar_tuples = new_read.cigartuples
                reference_middle = new_read.reference_start + round((
                    new_read.reference_length)/2)
                relative_pos = "left"

                if qseqs[new_read.qname]['poi'] < reference_middle:
                    # this is to fix the cigar library not handling correctly left mask if a left S mask already exists in the cigar string
                    # i.e. extra bases for masking are 6 to 35S103M, the new softmask should be (6+35)S97M
                    if cigar_tuples[0][0] == 4:
                        soft_clip_length += cigar_tuples[0][1]
                    masked_cigar_string = cigar_string.mask_left(
                        soft_clip_length)
                else:
                    masked_cigar_string = cigar_string.mask_right(
                        soft_clip_length)
                    relative_pos = "right"

                logger.opt(colors=True).info(
                    f" 🔅 <green>{new_read.qname}</green> has the dodgy base at {qseqs[new_read.qname]['poi']}, which is on the {relative_pos} of the mid point: {reference_middle}")

                update_message = f"▶️ <e>{cigar_string}</e> <r>➡️</r> <y>{masked_cigar_string}</y> | {qseqs[read.qname]['poi'] -  read.reference_start} | ▶️ query_length {new_read.query_length}"

                #Update new soft clipped cigar
                if len(masked_cigar_string) == new_read.query_length:
                    new_read.cigarstring = masked_cigar_string.cigar

                    #Adjust starting mapped location
                    new_read.reference_start += qseqs[new_read.qname]['poi'] - \
                        new_read.reference_start

                    logger.opt(colors=True).info(
                        f"{relative_pos.title()} Masked {update_message}")

                    if (qseqs[new_read.qname]['poi'] not in list(n_fixed_reads.keys())):
                        n_fixed_reads[qseqs[new_read.qname]['poi']] = 1
                    else:
                        n_fixed_reads[qseqs[new_read.qname]['poi']
                                      ] = n_fixed_reads[qseqs[new_read.qname]['poi']] + 1
                else:
                    logger.critical(f"❗ Excluded {update_message}")
        fixed_bam.write(new_read)
    bam.close()
    fixed_bam.close()
    print("Done!")
    os.system(f"samtools sort -@ {threads} -o {output_name} tmp.bam")
    os.system(f"samtools index {output_name}")
    # Python 8.8 missing_ok added
    pathlib.Path("tmp.bam").unlink()
    finished_message = f"Number of reads have been processed: {sum([int(v) for k,v in n_fixed_reads.items()])} {[f'{k}({v})' for k,v in n_fixed_reads.items()]}\nTotal reads found to be processed: {len(qseqs)}"
    # logger.info(finished_message)
    print(finished_message)


if __name__ == "__main__":
    logger.disable("__main__")
    soft_mask_dodgy_bases()
