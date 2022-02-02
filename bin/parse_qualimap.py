#!/usr/bin/env python3
import re
import subprocess
import csv
import sys

def parse_qualimap(input_file, sample_name):
    regexes = {
        'sequenced_bases' : r"number of sequenced bases\s*=\s*([\d,]+)",
        'reads': r"number of reads\s*=\s*([\d,]+)",
        'mapped_reads' : r"number of mapped reads\s*=\s*([\d,]+)",
        'coverage' : r"mean coverageData\s*=\s*([\d,\.X]+)"
    }
    _tsv_record = {'sample_name': sample_name}
    with open(input_file, 'r') as fh:
        _txt = fh.read()
        for k, v in regexes.items():
            _rs = re.findall(v,_txt)[0].replace(',','')
            _tsv_record.update({k : _rs})
        print(_tsv_record)
    
    output_file = f"{sample_name}.tsv"
    
    with open(output_file, "w") as fh:
        tsv = csv.DictWriter(fh, fieldnames=_tsv_record.keys())
        tsv.writeheader()
        tsv.writerow(_tsv_record)

def main():
    qualimap_file = sys.argv[1]
    sample_name = sys.argv[2]
    parse_qualimap(qualimap_file, sample_name)

if __name__ == "__main__":
    main()