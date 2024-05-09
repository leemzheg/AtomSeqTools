#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import pysam

def argparse_line():
    parser = argparse.ArgumentParser(description='Filter STAR alignments',
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', metavar='BAM_FILE', 
        help='the input bam file', required=True)
    required.add_argument('-o', metavar='BAM_FILE', 
        help='the output bam file', required=True)
    required.add_argument('--softclip-length', metavar='INT',
        help="filter alignment which 5' softclip length\n"
        "longer than {softclip-length} in read [Default:2]",
        required=True, type=int)
    required.add_argument('--paired', choices=['True', 'False'],
        help='paired input BAM', required=True)
    argv = vars(parser.parse_args())
    return argv

def pair_filter(input_bam_file, softclip_length, output_bam_file):
    with pysam.AlignmentFile(input_bam_file, "rb") as handle:
        with pysam.AlignmentFile(output_bam_file, "wb", 
            template=handle) as output:
            for read in handle:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.is_read2 and not read.is_supplementary:
                    continue
                ################################################################
                if read.is_read1:
                    if read.is_reverse:
                        if read.cigartuples[-1][0] == 4:
                            if read.cigartuples[-1][1] > softclip_length:
                                continue
                        output.write(read)
                    else:
                        if read.cigartuples[0][0] == 4:
                            if read.cigartuples[0][1] > softclip_length:
                                continue
                        output.write(read)
                else:
                    if read.is_reverse:
                        if read.cigartuples[0][0] == 4:
                            if read.cigartuples[0][1] > softclip_length:
                                continue
                        output.write(read)
                    else:
                        if read.cigartuples[-1][0] == 4:
                            if read.cigartuples[-1][1] > softclip_length:
                                continue
                        output.write(read)

def merge_filter(input_bam_file, softclip_length, output_bam_file):
    with pysam.AlignmentFile(input_bam_file, "rb") as handle:
        with pysam.AlignmentFile(output_bam_file, "wb", 
            template=handle) as output:
            for read in handle:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.is_reverse:
                    if read.cigartuples[-1][0] == 4:
                        if read.cigartuples[-1][1] > softclip_length:
                            continue
                    output.write(read)
                else:
                    if read.cigartuples[0][0] == 4:
                        if read.cigartuples[0][1] > softclip_length:
                            continue
                    output.write(read)

def main():
    argv = argparse_line()
    if argv['paired'] == 'True':
        pair_filter(argv['i'], argv['softclip_length'], argv['o'])
    else:
        merge_filter(argv['i'], argv['softclip_length'], argv['o'])

if __name__ == '__main__':
    main()
