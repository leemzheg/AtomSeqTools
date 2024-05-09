#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import json

def argparse_line():
    parser = argparse.ArgumentParser(description='Calculate primer depth', 
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('--intersect', metavar='INTERSECT',
        help='input intersect file between reads\n'
        'mapped region and primer region', required=True) 
    required.add_argument('--primer-bed', metavar='BED', 
        help='input primer bed file', required=True)
    required.add_argument('--primer-depth', metavar='PRIMER_DEPTH', 
        help='output primer depth file', required=True)
    required.add_argument('--read-primer', metavar='READ_PRIMER', 
        help='output read and primer corresponding files', required=True)
    required.add_argument('--distance', metavar='INT',
        help='distance between read mapped posithon\n'
        'and primer position [default:2]',
        required=True, type=int)
    argv = vars(parser.parse_args())
    return argv

def read_primer_bed(primer_bed_file):
    primer_num = {}
    primer_info = {}
    with open(primer_bed_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            primer_num[newline[3]] = 0
            primer_info[newline[3]] = newline
    return primer_num, primer_info

def read_intersect(intersect_file, primer_bed_file, 
                   distance, primer_depth_file, read_primer_file):
    
    primer_num, primer_info = read_primer_bed(primer_bed_file)
    
    read_primer = {}
    with open(intersect_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            mapped_start = int(newline[1])
            mapped_end = int(newline[2])
            read_id = newline[3]
            primer_start = int(newline[13])
            primer_end = int(newline[14])
            primer_name = newline[15]
            mapped_strand = newline[5]
            if mapped_strand == '+':
                if abs(mapped_start - primer_start) <= distance:
                    if read_id not in read_primer:
                        read_primer[read_id] = {}
                        read_primer[read_id][primer_name] = ''
                    else:
                        read_primer[read_id][primer_name] = ''
            else:
                if abs(mapped_end-primer_end) <= distance:
                    if read_id not in read_primer:
                        read_primer[read_id] = {}
                        read_primer[read_id][primer_name] = ''
                    else:
                        read_primer[read_id][primer_name] = ''
    
    with open(read_primer_file, 'w') as output:
        for read in read_primer:
            if len(read_primer[read]) > 1:
                print('This read has two primer: {0} {1}'.format(read, 
                                                      read_primer[read]))
            else:
                for primer in read_primer[read]:
                    primer_num[primer] += 1
                    output.write('{0}\t{1}\n'.format(read, primer))

    with open(primer_depth_file, 'w') as output:
        output.write('#Chromosome\tStart\tEnd\tPrimer'
                     '\tLength\tStrand\tAccession\tDepth\n')
        for primer in primer_num:
            output.write('{0}\t{1}\n'.format('\t'.join(primer_info[primer]),
                                             primer_num[primer]))

def main():
    argv = argparse_line()
    read_intersect(argv['intersect'], argv['primer_bed'], argv['distance'],
                   argv['primer_depth'], argv['read_primer'])

if __name__ == '__main__':
    main()
