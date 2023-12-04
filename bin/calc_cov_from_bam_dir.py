#!/usr/bin/env python

import subprocess
import os
import re
import argparse

def calc_cov(input_dir, output_file):
    files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
    with open(output_file, 'w') as f:
        for i in files:
           # if re.match('.*bai', i) is None: This replaces .bam with an empty string
            bin = re.sub('\.bam', '', i)
            cov = subprocess.check_output(f"samtools depth -a {input_dir}/{i} | awk '{{sum+=$3}} END {{print sum/NR}}'", shell=True)
            cov = cov.decode("utf-8").replace("\n", "")
            f.write(bin + "," + str(cov) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate coverage from BAM files')
    parser.add_argument('-i', '--input', help='Input directory containing BAM files', required=True)
    parser.add_argument('-o', '--output', help='Output coverage report file', required=True)
    args = parser.parse_args()

    calc_cov(args.input, args.output)
