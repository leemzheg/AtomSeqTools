#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import re
import gzip
import pysam

def argparse_line():
    parser = argparse.ArgumentParser(description='Extract pair reads '
        'from bam file to fastq file',
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('-bam', metavar='BAM', 
        help='the input bam file', required=True)
    required.add_argument('-fq1', metavar='FASTQ', 
        help='the output read1 file', required=True)
    required.add_argument('-fq2', metavar='FASTQ', 
        help='the output read2 file', required=True)
    argv = vars(parser.parse_args())
    return argv

def extract(bam_file, fq1_file, fq2_file):
    read1_sequence = {}
    read1_quality = {}
    read2_sequence = {}
    read2_quality = {}
    check_pair = {}
    with pysam.AlignmentFile(bam_file, "rb") as handle:
        for read in handle:
            quality = ''.join(map(lambda x: chr(x+33), 
                read.get_forward_qualities()))
            if read.qname not in check_pair:
                check_pair[read.qname] = {}
            if read.is_read1:
                check_pair[read.qname]['read1'] = ''
                read1_sequence[read.qname] = read.get_forward_sequence()
                read1_quality[read.qname] = quality
            elif read.is_read2:
                check_pair[read.qname]['read2'] = ''
                read2_sequence[read.qname] = read.get_forward_sequence()
                read2_quality[read.qname] = quality

    with gzip.open(fq1_file, 'wt') as read1:
        with gzip.open(fq2_file, 'wt') as read2:
            for read in check_pair:
                if len(check_pair[read]) == 2:
                    read1.write('@{0}\n{1}\n+\n{2}\n'.format(read, 
                                             read1_sequence[read],
                                             read1_quality[read]))
                    read2.write('@{0}\n{1}\n+\n{2}\n'.format(read, 
                                             read2_sequence[read],
                                             read2_quality[read]))

def main():
    argv = argparse_line()
    extract(argv['bam'], argv['fq1'], argv['fq2'])
if __name__ == '__main__':
    main()
