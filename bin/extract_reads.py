#!/usr/bin/env python3

import dnaio
import click
import arrow
from pathlib import Path


def is_included(start_time, start='2020-09-24T21:20:28Z', duration=30):
    _start_time = start_time.replace('start_time=','')
    _start_time = arrow.get(_start_time)
    time_start = arrow.get(start)
    time_end = time_start.shift(hours=duration)
    included = _start_time.is_between(time_start,time_end)
    return included

@click.command()
@click.argument('fastq-folder', required=False)
@click.option('--output-name', type=str)
@click.option('--start', default='2020-09-24T21:20:28Z', help='Start time of the sequencing run')
@click.option('--duration', type=int, default=30, help = 'Duration of the run (hours)')
def extract(fastq_folder, output_name, start, duration):
    fastq_files = [str(fq) for fq in Path(fastq_folder).iterdir(
    ) if fq.name.endswith('.fastq') or fq.name.endswith('.fastq.gz')]
    
    if fastq_folder is None:
        fastq_file = '/beegfs/Covid-19_Seq/cog13.coronahit/result.coronahit.2020926.cog13.guppy422.f55.rerun/articNcovNanopore_sequenceAnalysisMedaka_articGuppyPlex/NORW-COG13.coronaHiT_barcode01.fastq.gz'
    if len(fastq_files) >= 1:
        if not output_name:
            output_fastq = f"{Path(fastq_folder).name}.fastq.gz"
        else:
            output_fastq = f'{output_name}.fastq.gz'
        with dnaio.FastqWriter(output_fastq) as fo:
            for fastq_file in fastq_files:
                with dnaio.open(fastq_file) as fi:
                    for record in fi:
                        # seq_name,_,_,_,_,start_time,_,_ = record.name.split()
                        header = record.name.split()
                        seq_name = header[0]
                        start_time = header[5]
                        included_read = is_included(start_time, start = start, duration = duration)
                        if included_read:
                            fo.write(record)
                        else:
                            print(f'{seq_name}\t{start_time}\t{included_read}')

if __name__ == "__main__":
    extract()
