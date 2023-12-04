#!/usr/bin/env python

import argparse
import os
import shutil

def cov_filtered(input_file, input_dir, output_dir):
    print("Filter_coverage")
    res = []
    try:
        with open(input_file, 'r') as f:
            lines = [line.rstrip('\n') for line in f]
            for line in lines:
                coverage = int(line.split(',')[1].split(".")[0])
                if coverage > 30:
                    bin_name = line.split(',')[0]
                    fasta_file = os.path.join(input_dir, "{}.fa".format(bin_name))
                    output_fasta = os.path.join(output_dir, "{}.fa".format(bin_name))
                    res.append((fasta_file, output_fasta))
    except Exception as e:
        print("Error:", str(e))
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter fasta files based on coverage')
    parser.add_argument('-i', '--input', help='Input list file with coverage information', required=True)
    parser.add_argument('-d', '--input_dir', help='Input directory containing the fasta files', required=True)
    parser.add_argument('-o', '--output', help='Output directory to copy selected fasta files', required=True)
    args = parser.parse_args()

    result = cov_filtered(args.input, args.input_dir, args.output)

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Copy selected fasta files to the output directory
    for fasta_file, output_fasta in result:
        try:
            shutil.copy(fasta_file, output_fasta)
        except Exception as e:
            print("Error copying {}: {}".format(fasta_file, str(e)))
