#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import re
import json
from natsort import natsorted

def argparse_line():
    parser = argparse.ArgumentParser(description='Reads mapped summary', 
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('--fastp-json', metavar='JSON',
        help='fastp output json file', required=True)
    required.add_argument('--total-star-log', metavar='STAR_LOG',
        help='first star output log file', required=True)
    required.add_argument('--unique-star-log', metavar='STAR_LOG',
        help='second star output log file', required=True)
    required.add_argument('--total-primer-depth',
        metavar='TOTAL_PRIMER_DEPTH',
        help='total primer depth file', required=True)
    required.add_argument('--unique-primer-depth',
        metavar='UNIQUE_PRIMER_DEPTH',
        help='unique primer depth file', required=True)
    required.add_argument('--summary', metavar='SUMMARY', 
        help='the output summary file', required=True)
    argv = vars(parser.parse_args())
    return argv

def read_fastp_json(json_file):
    fastp = {}
    sample_name = json_file.split('/')[0]
    with open(json_file, 'r') as handle:
        js = json.load(handle)
        fastp[sample_name] = {}
        fastp[sample_name]['raw_reads'] = \
            js['summary']['before_filtering']['total_reads']
        fastp[sample_name]['clean_reads'] = \
            js['summary']['after_filtering']['total_reads']
        fastp[sample_name]['raw_reads_q30'] = \
            js['summary']['before_filtering']['q30_rate']
    
    return fastp

def read_star_log(star_log_file):
    star = {}
    sample_name = star_log_file.split('/')[0]
    with open(star_log_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if 'Number of input reads' in line:
                total_reads = int(newline[1])
            elif 'Uniquely mapped reads number' in line:
                unique_mapped_reads = int(newline[1])
            elif 'Number of reads mapped to multiple loci' in line:
                multi_loci_reads = int(newline[1])
            elif 'Number of reads mapped to too many loci' in line:
                many_loci_reads = int(newline[1])
            elif 'Number of chimeric reads' in line:
                chimeric_reads = int(newline[1])
        star[sample_name] = {}
        star[sample_name]['total_reads'] = total_reads
        star[sample_name]['mapped_reads'] = unique_mapped_reads + \
                multi_loci_reads + many_loci_reads + chimeric_reads
    
    return star

def calculate_uniformity(depth):
    mean = sum(depth) / len(depth)
    uniformity_threshold = mean * 0.2
    valid_uniformity_count = 0
    for i in depth:
        if i >= uniformity_threshold:
            valid_uniformity_count += 1
    uniformity = float(valid_uniformity_count) / len(depth)

    return sum(depth), mean, uniformity

def read_primer_depth(primer_depth_file):
    control = ['CHMP2A', 'GPI', 'RAB7A', 'VCP']
    primer_depth = {}
    sample_name = primer_depth_file.split('/')[0]
    all_depth = []
    ctrl_depth = []
    with open(primer_depth_file) as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if re.match('#', line):
                continue
            primer = newline[3]
            depth = int(newline[-1])
            gene_name = primer.split('_')[0]
            all_depth.append(depth)
            if gene_name in control:
                ctrl_depth.append(depth)
        depth_total, depth_mean, depth_uniformity = \
            calculate_uniformity(all_depth)
        if ctrl_depth:
            ctrl_depth_total, ctrl_depth_mean, ctrl_depth_uniformity = \
                calculate_uniformity(ctrl_depth)
        ############################################################
        primer_depth[sample_name] = {}
        primer_depth[sample_name]['depth_total'] = \
            int(depth_total)
        primer_depth[sample_name]['depth_mean'] = \
            int(depth_mean)
        primer_depth[sample_name]['depth_uniformity'] = \
            round(depth_uniformity*100, 2)
        if ctrl_depth:
            primer_depth[sample_name]['ctrl_depth_mean'] = \
                int(ctrl_depth_mean)
        else:
            primer_depth[sample_name]['ctrl_depth_mean'] = 'NA'
    
    return primer_depth

def summary(fastp_json_file, total_star_log_file, unique_star_log_file,
            total_primer_depth_file, unique_primer_depth_file,
            summary_file):
    
    fastp = read_fastp_json(fastp_json_file)
    total_star = read_star_log(total_star_log_file)
    unique_star = read_star_log(unique_star_log_file)
    total_primer_depth = read_primer_depth(total_primer_depth_file)
    unique_primer_depth = read_primer_depth(unique_primer_depth_file)

    with open(summary_file, 'w') as total_output:
        title = ('#Sample\tData_size(G)\tClean_rate%(PE)\t'
                 'Q30_rate%\tMapping_rate%(PE)\tUnique_read(PE)\t'
                 'Unique_percent%\tUnique_on_target_rate%\t'
                 'On_target_rate%\tUnique_mean_depth\t'
                 'Mean_depth\tUnique_uniformity%\tUniformity%\t'
                 'Unique_control_mean_depth\n')
        total_output.write(title)
        for sample in natsorted(fastp.keys()):
            data_size = round(
                    int(fastp[sample]['raw_reads'])*150/1000000000, 2)
            raw_read = int(fastp[sample]['raw_reads']/2)
            clean_read = int(fastp[sample]['clean_reads']/2)
            clean_rate = round(clean_read/raw_read*100, 2)
            q30_rate = round(fastp[sample]['raw_reads_q30']*100, 2)
            mapping_rate = round(float(total_star[sample]['mapped_reads'])/ \
                           (total_star[sample]['total_reads'])*100, 2)
            unique_read = int(unique_star[sample]['total_reads'])
            unique_percent = round(unique_read/clean_read*100, 2)
            unique_on_target_rate = round(
                unique_primer_depth[sample]['depth_total']/ \
                unique_star[sample]['mapped_reads']*100, 2)
            on_target_rate = round(
                total_primer_depth[sample]['depth_total']/ \
                total_star[sample]['mapped_reads']*100, 2)
            unique_mean_depth = unique_primer_depth[sample]['depth_mean']
            mean_depth = total_primer_depth[sample]['depth_mean']
            unique_uniformity = \
                    unique_primer_depth[sample]['depth_uniformity']
            uniformity = total_primer_depth[sample]['depth_uniformity']
            unique_control_mean_depth = \
                    unique_primer_depth[sample]['ctrl_depth_mean']
            summary = (f'{sample}\t{data_size}\t{clean_rate}\t'
                       f'{q30_rate}\t{mapping_rate}\t{unique_read}\t'
                       f'{unique_percent}\t{unique_on_target_rate}\t'
                       f'{on_target_rate}\t{unique_mean_depth}\t'
                       f'{mean_depth}\t{unique_uniformity}\t{uniformity}\t'
                       f'{unique_control_mean_depth}\n')
            total_output.write(summary)

def main():
    argv = argparse_line()
    summary(argv['fastp_json'],
            argv['total_star_log'], argv['unique_star_log'],
            argv['total_primer_depth'], argv['unique_primer_depth'],
            argv['summary'])
if __name__ == '__main__':
    main()
