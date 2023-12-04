#!/usr/bin/env python

import subprocess
import argparse
import os
import shutil
import re

def get_best(input_file, output_dir, bam_dir):
    res = []
    print("Filter by MAG quality")
    test = """awk '$2>=80 && $3<=20 {{ print $1}}' {} """.format(input_file)
    proc = subprocess.Popen(test, shell=True, stdout=subprocess.PIPE)
    stdout_value = proc.communicate()[0]
    stdout = re.sub('\'', '', re.sub('b', '', str(stdout_value), 1))
    for item in stdout.split('\\n'):
        if item != '':
            bam_file = os.path.join(bam_dir, "{}.bam".format(item))
            bai_file = os.path.join(bam_dir, "{}.bam.bai".format(item))
            res.append((bam_file, bai_file))
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select best quality BAM files')
    parser.add_argument('-i', '--input', help='Input quality report TSV file', required=True)
    parser.add_argument('-o', '--output', help='Output directory', required=True)
    parser.add_argument('-b', '--bamdir', help='Directory containing the original BAM files', required=True)
    args = parser.parse_args()
    
    # create output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
	
    result = get_best(args.input, args.output, args.bamdir)
  #  print(result)

# Copy selected BAM files to the output directory
    for bam, bai in result:
        shutil.copy(bam, args.output)
        shutil.copy(bai, args.output)
