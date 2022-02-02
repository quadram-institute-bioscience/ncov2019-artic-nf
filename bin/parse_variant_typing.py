#!/usr/bin/env python3

import pandas as pd
import sidetable
import click
import pathlib
from discord_webhook import DiscordWebhook as DWH

test_url = ""
qib_discord_url = ""

discord_urls = [test_url, qib_discord_url]

@click.command()
@click.option('--name', default = "pangolin")
@click.argument('report-file', required=False)
def parsing(report_file,name):
    # if report_file is None:
    #     report_file = "/share/civet_report/20201216.B.1.1.7/pangolin/lineage_report.csv"
    assert pathlib.Path(report_file).exists(), f"{report_file} is not existed"
    report = pd.read_csv(report_file)
    report = report[["sampleID", "type"]]
    lineage_report = report.stb.freq(['type'])[['type', 'count']]
    lineage_report.columns = ['type', 'no of samples']
    lineage_report = lineage_report.round(2)
    username = "nCovid19-pipeline"
    wh = DWH(discord_urls, username = username)
    wh.content = f":point_right: **Variants of interest (by VCF typing) found for {name}** \n```{str(lineage_report)}\n```"
    wh.execute()
    
if __name__ == "__main__":
    parsing()
