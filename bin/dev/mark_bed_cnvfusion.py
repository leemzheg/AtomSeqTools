#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap

def argparse_line():
	parser = argparse.ArgumentParser(description='', 
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--input', metavar='', 
		help='', required=True)
	parser.add_argument('--cnv-gene', metavar='', 
		help='', required=True)
	parser.add_argument('--output', metavar='', 
		help='', required=True)
	argv = vars(parser.parse_args())
	return argv

def mark(inputfile, cnv_gene_file, outputfile):
	cnv_gene = []
	for line in open(cnv_gene_file, 'r'):
		cnv_gene.append(line.strip())

	with open(outputfile, 'w') as output:
		for line in open(inputfile, 'r'):
			newline = line.strip().split('\t')
			geneset = newline[3].split(',')
			newgeneset = []
			if '_' in newline[0]:
				continue
			for gene in geneset:
				if gene in cnv_gene:
					newgeneset.append(gene)
			if len(newgeneset) > 1:
				print(newline, newgeneset)
				exit()
			if len(newgeneset) != 0:
				newline[3] = ''.join(newgeneset) + ':CNV'
				output.write('{0}\n'.format('\t'.join(newline)))
			if len(newgeneset) == 0:
				output.write('{0}\n'.format('\t'.join(newline)))

def main():
	argv = argparse_line()
	mark(argv['input'], argv['cnv_gene'], argv['output'])

if __name__ == '__main__':
	main()
