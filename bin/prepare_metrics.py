#!/usr/bin/env python3
import re
import pandas as pd
import sys

def parse_qualimap(mapping_metric, count_snps, qc_metric, run_name):
    mapping = pd.read_csv(mapping_metric)
    mapping.rename(columns=lambda x: x.strip(), inplace=True)
    count_snps = pd.read_csv(count_snps)
    count_snps.rename(columns=lambda x: x.strip(), inplace=True)
    qc = pd.read_csv(qc_metric)
    qc.rename(columns=lambda x: x.strip(), inplace=True)
    
    all_metrics = [df.set_index(['sample_name']) for df in [mapping, count_snps, qc]]
    all_metrics = pd.concat(all_metrics, axis=1).reset_index()
    all_metrics.rename(columns={'index':'sample_name'}, inplace=True)

    # all_metrics['central_sample'] = 'NORW-' + all_metrics['sample_name'].str.split("_").str[0]
    all_metrics['central_sample'] = all_metrics['sample_name'].str.split("_").str[0]
    all_metrics['hi_qc_pass'] = all_metrics['pct_N_bases'].apply(
        lambda x: 'True' if x < 10.0 and x > 0 else 'False')
    all_metrics.drop(['fasta', 'bam'], axis=1, inplace=True)
    central_sample = all_metrics.pop('central_sample')
    all_metrics.insert(0, "central_sample",central_sample)
    all_metrics.to_csv(f"{run_name}_metrics.csv", index = False)

def main():
    mapping_metric = sys.argv[1]
    count_snps = sys.argv[2]
    qc_metric = sys.argv[3]
    run_name = sys.argv[4]
    parse_qualimap(mapping_metric, count_snps, qc_metric, run_name)

if __name__ == "__main__":
    main()
