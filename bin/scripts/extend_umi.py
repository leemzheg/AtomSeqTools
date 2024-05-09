#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def argparse_line():
    parser = argparse.ArgumentParser(description=
        "Extract part of the read2 5' seq as UMI extend to read ID", 
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', metavar='FASTQ', 
        help='read1 input file', required=True)
    required.add_argument('-o', metavar='FASTQ', 
        help='read1 output file', required=True)
    required.add_argument('-I', metavar='FASTQ', 
        help='read2 input file', required=True)
    required.add_argument('-O', metavar='FASTQ', 
        help='read2 output file', required=True)
    required.add_argument('--extend-length', metavar='INT', 
        help="specify UMI extend length [Default:11]", required=True, type=int)
    argv = vars(parser.parse_args())
    return argv

def extend(r1_input_file, r1_output_file, r2_input_file, r2_output_file, 
            extend_length):
    with gzip.open(r1_output_file, 'wt', compresslevel=4) as r1_output:
        with gzip.open(r2_output_file, 'wt', compresslevel=4) as r2_output:
            with gzip.open(r1_input_file, 'rt') as handle1:
                with gzip.open(r2_input_file, 'rt') as handle2:
                    for (title1, seq1, qual1), (title2, seq2, qual2) in \
                        zip(FastqGeneralIterator(handle1), \
                        FastqGeneralIterator(handle2)):
                        extend_umi = seq2[:extend_length]
                        title1 = title1.split()[0] + extend_umi
                        title2 = title2.split()[0] + extend_umi
                        r1_output.write('@{0}\n{1}\n+\n{2}\n'.format(
                            title1, seq1, qual1))
                        r2_output.write('@{0}\n{1}\n+\n{2}\n'.format(
                            title2, seq2, qual2))

def main():
    argv = argparse_line()
    extend(argv['i'], argv['o'], argv['I'], argv['O'], argv['extend_length']) 
if __name__ == '__main__':
    main()
