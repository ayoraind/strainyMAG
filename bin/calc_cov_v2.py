#!/usr/bin/env python

import subprocess
import re
import argparse

# usage: calc_cov_v2.py -i bam_file -o output_file

def calc_cov(bam_file, output_file):
    bin = re.sub(r'\.bam', '', bam_file)
    cov = subprocess.check_output(
        f"samtools depth -a {bam_file} | awk '{{sum+=$3}} END {{print sum/NR}}'", shell=True)
    cov = cov.decode("utf-8").replace("\n", "")
    with open(output_file, 'w') as output:
        output.write(f"{bin},{cov}\n")

def main():
    parser = argparse.ArgumentParser(description='Calculate coverage from a BAM file.')
    parser.add_argument('-i', '--input', help='Input BAM file path', required=True)
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    args = parser.parse_args()

    calc_cov(args.input, args.output)

if __name__ == "__main__":
    main()
