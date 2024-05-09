#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse, textwrap

parser = argparse.ArgumentParser(
    description="", formatter_class=argparse.RawTextHelpFormatter
)
g = parser.add_argument_group("required arguments")
g.add_argument("-primer", metavar="i", help="i", required=True)
g.add_argument("-outdir", metavar="o", help="o", required=True)
argv = parser.parse_args()

primerbed = argv.primer
outd = argv.outdir
with open(primerbed) as f, open(f"{outd}/zz.tmp.xls", "w") as fi:
    for line in f:
        ls = line.strip().split("\t")
        if ls[5] == "+":
            fi.write("{}\t{}\t{}\n".format(ls[0], ls[2], str(int(ls[2]) + 80)))
        elif ls[5] == "-":
            fi.write("{}\t{}\t{}\n".format(ls[0], str(int(ls[1]) - 80), ls[1]))
        else:
            print("Wrong: BED format has wrong, please check")
            exit()
