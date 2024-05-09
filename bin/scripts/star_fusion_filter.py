#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import re
import os
import ast
import json
import pysam
import gffutils

def argparse_line():
    parser = argparse.ArgumentParser(description='Calculate fusion frequency '
        'and fusion per primer frequency', 
        formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('required arguments')
    required.add_argument('--prediction', metavar='PREDICTION', 
        help='star-fusion.fusion_predictions.tsv', 
        required=True)
    required.add_argument('--coding-effect', metavar='CODING_EFFECT', 
        help='star-fusion.fusion_predictions.abridged.coding_effect.tsv', 
        required=True)
    required.add_argument('--ctat-cancer-introns', metavar='CANCER_INTRONS', 
        help='*.cancer.introns', 
        required=True)
    required.add_argument('--ctat-introns', metavar='INTRONS', 
        help='*.introns', 
        required=True)
    required.add_argument('--primer-depth', metavar='PRIMER_DEPTH', 
        help='primer_depth.txt', 
        required=True)
    required.add_argument('--read-primer', metavar='INTERSECT',
        help='read_primer.txt', 
        required=True)
    required.add_argument('--genome-lib-dir', metavar='PATH',
        help='directory containing genome lib', 
        required=True)
    required.add_argument('--input-bam', metavar='BAM',
        help='Aligned.sortedByCoord.out.bam', 
        required=True)
    required.add_argument('--output-bam', metavar='BAM',
        help='*.intarget_junction_reads.bam', 
        required=True)
    required.add_argument('--output-bedpe', metavar='BEDPE',
        help='*.fusion_frequency.bedpe', 
        required=True)
    required.add_argument('--igv-json', metavar='JSON',
        help='*.fusion_frequency.trackConfigs.json', 
        required=True)
    required.add_argument('--primer-frequency', metavar='PRIMER_PREQUENCY', 
        help='star-fusion.primer_frequency.xls', 
        required=True)
    required.add_argument('--fusion-frequency', metavar='FUSION_PREQUENCY', 
        help='star-fusion.fusion_frequency.xls', 
        required=True)
    required.add_argument('--min-junction-reads', metavar='MIN_JUNCTION_READS', 
        help='minimum number of junction reads required [default:5]', type=int, 
        required=True)
    required.add_argument('--min-fusion-frequency', metavar='MIN_FREQUENCY', 
        help='minimum fusion frequency [default:0.01]', type=float,
        required=True)
    required.add_argument('--min-primer-depth', metavar='MIN_PRIMER_DEPTH', 
        help='minimum fusion all primer depth [default:50]', type=int,
        required=True)
    required.add_argument('--fusion-filter', metavar='FUSION_FILTER', 
        help='star-fusion.fusion_filter.xls', 
        required=True)
    argv = vars(parser.parse_args())
    return argv

def read_primer_depth(primer_depth_file):
    depth = {}
    kit_accession = []
    with open(primer_depth_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if re.match('#', line):
                continue
            for accession in newline[6].split(','):
                if accession not in kit_accession and accession != 'NA':
                    kit_accession.append(accession)
            primer_name = newline[3]
            primer_depth = newline[-1]
            depth[primer_name] = int(primer_depth)
    return depth, kit_accession

def read_read_primer(read_primer_file):
    read_primer = {}
    with open(read_primer_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            read = newline[0]
            primer_name = newline[1]
            read_primer[read] = primer_name
    return read_primer

def read_coding_effect(coding_effect_file):
    prot_fusion_type = {}
    with open(coding_effect_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if re.match('#', line):
                continue
            transcript = newline[0] + '_' + newline[7] + '_' + newline[9]
            prot_fusion_type[transcript] = newline[21]
    return prot_fusion_type

def read_oncokb_gene():
    oncokb_accession = []
    tsv_dir = os.path.dirname(os.path.abspath(__file__))
    with open('{0}/cancerGeneList.tsv'.format(tsv_dir), 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if re.match('Hugo', line):
                continue
            if newline[5] and newline[5] != 'null':
                accession = newline[5].split('.')[0]
                oncokb_accession.append(accession)
    return oncokb_accession

def calculate_primer_frequency(junction_reads, read_primer):
    fusion_primer_depth = {}
    fusion_primer_read = {}
    fusion_intarget_read = []
    fusion_outoftarget_read = []
    for read in junction_reads.split(','):
        read_name = read.split('@')[-1]
        if read_name in read_primer:
            fusion_intarget_read.append(read_name)
            primer = read_primer[read_name]
            if primer not in fusion_primer_depth:
                fusion_primer_depth[primer] = 1
                fusion_primer_read[primer] = [read_name] 
            else:
                fusion_primer_depth[primer] += 1
                fusion_primer_read[primer].append(read_name)
        else:
            fusion_outoftarget_read.append(read_name)
    return fusion_primer_depth, fusion_primer_read, \
           fusion_intarget_read, fusion_outoftarget_read

def find_exon_number(break_point, genome_lib_dir, 
                     kit_accession, oncokb_accession):

    db = gffutils.FeatureDB('{0}/ref_annot.gtf.sqldb'.format(genome_lib_dir), 
         keep_order=True)

    break_point = break_point.split(':')
    region = '{0}:{1}-{1}'.format(break_point[0], break_point[1])
    strand = break_point[2]
    
    candidate_transcript = {'Kit_Select':'', 'OncoKB_Select':'',
                            'MANE_Select':'', 'RefSeq_Select':''}
    for i in db.region(region, strand=strand, featuretype=['exon', 'intron']):
        transcript = ''.join(i.attributes['transcript_id']).split('.')[0]
        if i.featuretype == 'exon':
            number = ''.join(i.attributes['exon_number'])
        else:
            number = i.attributes['exon_number'][0]
        #####################################################################
        if transcript in kit_accession:
            candidate_transcript['Kit_Select'] = \
                transcript + ':' + i.featuretype + number
        if transcript in oncokb_accession:
            candidate_transcript['OncoKB_Select'] = \
                transcript + ':' + i.featuretype + number
        if 'tag' in i.attributes:
            if 'MANE Select' in i.attributes['tag']:
                candidate_transcript['MANE_Select'] = \
                    transcript + ':' + i.featuretype + number
            if 'RefSeq Select' in i.attributes['tag']:
                candidate_transcript['RefSeq_Select'] = \
                    transcript + ':' + i.featuretype + number 
                    
    exon_number = []
    for select in candidate_transcript:
        if candidate_transcript[select]:
            exon_number.append(select + ':' + candidate_transcript[select])
    
    if not exon_number:
        primary_exon_num = 'NA'
        exon_number = 'NA'
    else:
        primary_exon_num = exon_number[0].split(':')[-1]
        exon_number = ','.join(exon_number)
    
    return primary_exon_num, exon_number

def read_ctat_splicing(ctat_cancer_introns, ctat_introns):
    target = {'METx14del':'MET_Exon14_Skipping',
              'EGFRvIII':'EGFR_Exon2-7_Skipping'}

    exon_skipping = {}
    exon_range = {}
    with open(ctat_cancer_introns, 'r') as handle:
        for line in handle:
            newline = line.rstrip().split('\t')
            if re.match('intron', line):
                continue
            variant_name = newline[7]
            read_num = int(newline[3])
            loc = newline[0]
            chromosome = loc.split(':')[0]
            start = int(loc.split(':')[-1].split('-')[0])
            end = int(loc.split(':')[-1].split('-')[-1])
            if variant_name in target:
                exon_skipping[target[variant_name]] = {}
                exon_skipping[target[variant_name]]['readnum'] = read_num
                exon_skipping[target[variant_name]]['loc'] = loc
                exon_skipping[target[variant_name]]['totalread'] = 0
                exon_range[target[variant_name]] = {}
                exon_range[target[variant_name]][chromosome] = [start, end]
   
    with open(ctat_introns, 'r') as handle:
        for line in handle:
            newline = line.rstrip().split('\t')
            if re.match('intron', line):
                continue
            chromosome = newline[0].split(':')[0]
            start = int(newline[0].split(':')[-1].split('-')[0])
            end = int(newline[0].split(':')[-1].split('-')[-1])
            uniq_mapped = int(newline[3])
            for skipping in exon_range:
                if chromosome in exon_range[skipping]:
                    if exon_range[skipping][chromosome][0] <= start and \
                        exon_range[skipping][chromosome][1] >= end:
                        exon_skipping[skipping]['totalread'] += uniq_mapped
    return exon_skipping

def extract_fusion_intarget_read_align(input_bam_file, output_bam_file, 
    all_fusion_intarget_read, genome_lib_dir, igv_json_file, skipping):

    all_read_id = []
    with pysam.AlignmentFile(input_bam_file, "rb") as handle:
        with pysam.AlignmentFile(output_bam_file, "wb",
            template=handle) as output:
            for gene in skipping:
                skipping_reads = handle.fetch(skipping[gene]['chr'], skipping[gene]['start'], skipping[gene]['end'])
                for read in skipping_reads:
                    if read not in all_read_id:
                        output.write(read)
                        all_read_id.append(read)
            for read in handle:
                if read.query_name in all_fusion_intarget_read:
                    if read not in all_read_id:
                        output.write(read)
                        all_read_id.append(read)

    sorted_bam = '{0}.sorted.bam'.format(output_bam_file.split('.bam')[0])
    pysam.sort("-o", sorted_bam, output_bam_file)
    pysam.index(sorted_bam)
    
    track_config = [
        {"name": "Refseq Genes",
         "format": "refgene",
         "url": "{0}/ncbiRefSeq.sorted.txt".format(genome_lib_dir)
        },
        {
        "name": "Alignments",
        "url": "{0}".format(sorted_bam),
        "indexURL": "{0}.bai".format(sorted_bam),
        "showSoftClips": "true"
        }
    ]
    
    with open(igv_json_file, 'w') as handle:
        json.dump(track_config, handle)

def prediction_filter(input_file, coding_effect_file,
        ctat_cancer_introns, ctat_introns,
        primer_depth_file, read_primer_file, genome_lib_dir, 
        input_bam_file, output_bam_file, output_bedpe_file, igv_json_file, 
        primer_frequency_file, fusion_frequency_file, min_junction_reads, 
        min_fusion_frequency, min_primer_depth, fusion_filter_file):
    
    primer_depth, kit_accession = read_primer_depth(primer_depth_file)
    read_primer = read_read_primer(read_primer_file)
    prot_fusion_type = read_coding_effect(coding_effect_file)
    exon_skipping = read_ctat_splicing(ctat_cancer_introns, ctat_introns)
    oncokb_accession = read_oncokb_gene()

    output_bedpe = open(output_bedpe_file, 'w')

    primer_frequency = open(primer_frequency_file, 'w')
    primer_frequency.write('#FusionName\tPrimer\t'
                           'InTargetJunctionReadCount\t'
                           'PrimerReadCount\t'
                           'Frequency\tLeftGene\tLeftBreakpoint\t'
                           'RightGene\tRightBreakpoint\t'
                           'PrimerJunctionRead\tAnnots\n')
    
    fusion_frequency = open(fusion_frequency_file, 'w')
    fusion_frequency.write('#FusionName\tJunctionReadCount\t'
                           'InTargetJunctionReadCount\t'
                           'OutOfTargetJunctionReadCount\t'
                           'AllPrimerReadCount\t'
                           'InTargetJunctionReadFrequency\t'
                           'ProtFusionType\t'
                           'LeftGene\tLeftBreakpoint\t'
                           'LeftExonNumber\t'
                           'RightGene\tRightBreakpoint\t'
                           'RightExonNumber\t'
                           'Annots\t'
                           'InTargetJunctionRead\t'
                           'OutOfTargetJunctionRead\n')

    fusion_filter = open(fusion_filter_file, 'w')
    fusion_filter.write('#FusionName\t'
                        'ReadNumber\t'
                        'PrimerReadNumber\t'
                        'FusionFrequency\t'
                        'ProtFusionType\t'
                        'LeftGene\tLeftBreakpoint\t'
                        'LeftExonNumber\t'
                        'RightGene\tRightBreakpoint\t'
                        'RightExonNumber\t'
                        'Annots\n')

    check_fusion_splicing = {}
    all_fusion_intarget_read = []
    with open(input_file, 'r') as handle:
        for line in handle:
            newline = line.strip().split('\t')
            if re.match('#', line):
                continue
            ###########################################################
            fusion_name = newline[0]
            if fusion_name not in check_fusion_splicing:
                check_fusion_splicing[fusion_name] = 0
            else:
                check_fusion_splicing[fusion_name] += 1
            ###########################################################
            junction_read_count = int(newline[1])
            left_gene = newline[6].split('^')[0]
            left_break_point = newline[7]
            right_gene = newline[8].split('^')[0]
            right_break_point = newline[9]
            junction_reads = newline[10]
            annots = newline[18]
            transcript = fusion_name + '_' + left_break_point + \
                         '_' + right_break_point
            ###########################################################
            fusion_primer_depth, fusion_primer_read, \
            fusion_intarget_read, fusion_outoftarget_read \
            = calculate_primer_frequency(junction_reads, read_primer)
            if not fusion_primer_depth:
                continue
            ######################primer frequency#####################
            fusion_all_primer_depth = 0
            for primer in fusion_primer_depth:
                fusion_all_primer_depth += primer_depth[primer]
                primer_fre = round(fusion_primer_depth[primer] / \
                             primer_depth[primer]*100, 4)
                primer_frequency.write(
                                '{0}_splicing{1}\t{2}\t{3}\t'
                                '{4}\t{5}\t{6}\t{7}\t{8}\t'
                                '{9}\t{10}\t{11}\n'.format(
                                fusion_name,
                                check_fusion_splicing[fusion_name],
                                primer,
                                fusion_primer_depth[primer],
                                primer_depth[primer],
                                primer_fre,
                                left_gene,
                                left_break_point,
                                right_gene,
                                right_break_point,
                                ','.join(fusion_primer_read[primer]),
                                annots))
            #########################################################
            left_primary_exon_number, left_exon_number = \
                               find_exon_number(left_break_point, 
                                                genome_lib_dir, 
                                                kit_accession,
                                                oncokb_accession)
            right_primary_exon_number, right_exon_number = \
                               find_exon_number(right_break_point, 
                                                genome_lib_dir, 
                                                kit_accession,
                                                oncokb_accession)
            #########################################################
            if left_primary_exon_number != 'NA' and \
                right_primary_exon_number != 'NA':
                fusion_name = '{0}({1})-{2}({3})'.format(left_gene,
                                                    left_primary_exon_number,
                                                    right_gene,
                                                    right_primary_exon_number)
            else:
                fusion_name = '{0}-{1}'.format(left_gene, right_gene)
            ##########################bed.pe##########################
            #left_chr = left_break_point.split(':')[0]
            #left_position = left_break_point.split(':')[1]
            #right_chr = right_break_point.split(':')[0]
            #right_position = right_break_point.split(':')[1]
            #output_bedpe.write('{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n'.format(
            #                    left_chr, left_position, right_chr, 
            #                    right_position, fusion_name))
            #####################fusion frequency#####################
            all_fusion_intarget_read += fusion_intarget_read
            fusion_intarget_fre = len(fusion_intarget_read) / \
                                  fusion_all_primer_depth
            fusion_frequency.write('{0}\t{1}\t{2}\t{3}\t'
                                   '{4}\t{5}\t{6}\t{7}\t'
                                   '{8}\t{9}\t{10}\t{11}\t'
                                   '{12}\t{13}\t{14}\t{15}\n'.format(
                                   fusion_name,
                                   junction_read_count,
                                   len(fusion_intarget_read),
                                   len(fusion_outoftarget_read),
                                   fusion_all_primer_depth,
                                   round(fusion_intarget_fre*100, 4),
                                   prot_fusion_type[transcript],
                                   left_gene,
                                   left_break_point,
                                   left_exon_number,
                                   right_gene,
                                   right_break_point,
                                   right_exon_number,
                                   annots,
                                   ','.join(fusion_intarget_read),
                                   ','.join(fusion_outoftarget_read)))
            #####################fusion filter#######################
            if len(fusion_intarget_read) < min_junction_reads:
                continue
            if fusion_intarget_fre <= min_fusion_frequency:
                continue
            if fusion_all_primer_depth <= min_primer_depth:
                continue
            if prot_fusion_type[transcript] == 'FRAMESHIFT':
                continue
            fusion_filter.write('{0}\t{1}\t{2}\t{3}\t{4}\t'
                                '{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.format(
                                fusion_name,
                                len(fusion_intarget_read),
                                fusion_all_primer_depth,
                                round(fusion_intarget_fre*100, 4),
                                prot_fusion_type[transcript],
                                left_gene,
                                left_break_point,
                                left_exon_number,
                                right_gene,
                                right_break_point,
                                right_exon_number,
                                annots))
            #########################bed.pe##########################
            left_chr = left_break_point.split(':')[0]
            left_position = left_break_point.split(':')[1]
            right_chr = right_break_point.split(':')[0]
            right_position = right_break_point.split(':')[1]
            output_bedpe.write('{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n'.format(
                                left_chr, left_position, right_chr, 
                                right_position, fusion_name))
    ########################################################################
    skipping = {}
    if exon_skipping:
        for i in exon_skipping:
            frequency = round(
                exon_skipping[i]['readnum']/
                exon_skipping[i]['totalread']*100, 2)
            left_break_point = exon_skipping[i]['loc'].split('-')[0]
            right_break_point = exon_skipping[i]['loc'].split(':')[0] + ':' +\
                exon_skipping[i]['loc'].split('-')[-1]
            if i == 'MET_Exon14_Skipping':
                fusion_filter.write('{0}\t{1}\t{2}\t{3}\tNA\t'
                                    'MET\t{4}\t'
                                    'Kit_Select:NM_000245:exon14\t'
                                    'MET\t{5}\t'
                                    'Kit_Select:NM_000245:exon14\t'
                                    'NA\n'.format(
                                    i, exon_skipping[i]['readnum'],
                                    exon_skipping[i]['totalread'],
                                    frequency,
                                    left_break_point,
                                    right_break_point
                                ))
                #########################bed.pe##########################
                left_chr = left_break_point.split(':')[0]
                left_position = left_break_point.split(':')[1]
                right_chr = right_break_point.split(':')[0]
                right_position = right_break_point.split(':')[1]
                output_bedpe.write('{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n'.format(
                                    left_chr, left_position, right_chr, 
                                    right_position, i))
                skipping['MET'] = {'chr':left_chr, 'start':int(left_position), 'end':int(right_position)}
            elif i == 'EGFR_Exon2-7_Skipping':
                fusion_filter.write('{0}\t{1}\t{2}\t{3}\tNA\t'
                                    'EGFR\t{4}\t'
                                    'Kit_Select:NM_005228:exon2-7\t'
                                    'EGFR\t{5}\t'
                                    'Kit_Select:NM_005228:exon2-7\t'
                                    'NA\n'.format(
                                    i, exon_skipping[i]['readnum'],
                                    exon_skipping[i]['totalread'],
                                    frequency,
                                    left_break_point,
                                    right_break_point
                                ))
                #########################bed.pe##########################
                left_chr = left_break_point.split(':')[0]
                left_position = left_break_point.split(':')[1]
                right_chr = right_break_point.split(':')[0]
                right_position = right_break_point.split(':')[1]
                output_bedpe.write('{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n'.format(
                                    left_chr, left_position, right_chr, 
                                    right_position, i))
                skipping['EGFR'] = {'chr':left_chr, 'start':int(left_position), 'end':int(right_position)}
    ########################################################################
    extract_fusion_intarget_read_align(input_bam_file, output_bam_file,
                                  all_fusion_intarget_read, genome_lib_dir,
                                  igv_json_file, skipping)
    ########################################################################
    primer_frequency.close()
    fusion_frequency.close()
    fusion_filter.close()
    output_bedpe.close()

def main():
    argv = argparse_line()
    prediction_filter(argv['prediction'], argv['coding_effect'],
                      argv['ctat_cancer_introns'], argv['ctat_introns'],
                      argv['primer_depth'], argv['read_primer'],
                      argv['genome_lib_dir'], argv['input_bam'],
                      argv['output_bam'], argv['output_bedpe'],
                      argv['igv_json'], argv['primer_frequency'],
                      argv['fusion_frequency'], argv['min_junction_reads'],
                      argv['min_fusion_frequency'], argv['min_primer_depth'],
                      argv['fusion_filter'])

if __name__ == '__main__':
    main()
