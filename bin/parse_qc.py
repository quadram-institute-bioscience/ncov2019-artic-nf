#!/usr/bin/env python3

import pandas as pd
import click
import pathlib
from discord_webhook import DiscordWebhook as DWH

test_url = ""
discord_urls = [test_url]

@click.command()
@click.argument('qc-file', required=False)
def parsing(qc_file):
    if qc_file is None:
        qc_file = "/beegfs/Covid-19_Seq/result.illumina.20201014/NORW-20201014.qc.csv"
    assert pathlib.Path(qc_file).exists(), f"{qc_file} is not existed"
    qc = pd.read_csv(qc_file)
    # qc.rename(columns = {'qc_pass':'COGUK_pass'}, inplace=True)
    # qc['COGUK_pass'].replace({True: 'Yes', False: 'No'}, inplace=True)
    qc['COGUK_pass'] = qc['pct_covered_bases'].apply(
        lambda x: 'Yes' if x >= 50 else 'No')
    qc['GISAID_pass'] = qc['pct_covered_bases'].apply(
        lambda x: 'Yes' if x >= 90 else 'No')
    qc_pass = pd.crosstab(qc.COGUK_pass, qc.GISAID_pass,
                      margins=True, margins_name="Total")
    low_full_reads = qc[qc.full_mapped_reads < 100][["sample_name","full_mapped_reads"]]
    low_full_reads.sort_values(by = ["sample_name"], inplace=True)
    # print(qc_pass)
    # print(low_full_reads.set_index("sample_name").to_dict())
    #QC
    wh = DWH(discord_urls)
    wh.content = f":point_right: **QC pass**\n```{str(qc_pass)}```"
    wh.execute()
    # Low fully mapped reads
    low_full_reads_list = []
    for index, row in low_full_reads.iterrows():
        low_full_reads_list.append(
            f'{row["sample_name"]} {row["full_mapped_reads"]:<5}')
    low_full_reads_context = "\n".join(low_full_reads_list)
    wh.content = f":point_right: **{len(low_full_reads_list)} samples with < 100 fully mapped reads**\n```{low_full_reads_context}```"
    wh.execute()

if __name__ == "__main__":
    parsing()
