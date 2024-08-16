#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/03/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os
import argparse, textwrap
import pandas as pd
from multiprocessing import Pool
from functools import partial
import datetime

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="Variant QC threshold software. From excel to excel",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="cfDNA recommand threshold:\n  depth >= 500\n  alt >= 5\n  frequency >= 0.5\n"
        "Tissue recommand threshold:\n  depth >= 300\n  alt >= 5\n  frequency >= 1\n",
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-snvIndels", metavar="STR", help="Input 10_snv_indels.xls", required=True
    )
    required.add_argument(
        "-depth",
        metavar="INT",
        help="Threshold about mutation site total depth (>=depth). [default: 300]",
        type=int,
        default=300,
    )
    required.add_argument(
        "-alt",
        metavar="INT",
        help="Threshold about mutation alternative allele (>=alt). [default: 5]",
        type=int,
        default=5,
    )
    required.add_argument(
        "-frequency",
        metavar="FLOAT",
        help="Threshold about mutation site frequency (>=frequency). [default: 1]",
        type=float,
        default=1,
    )
    required.add_argument("-out", metavar="STR", help="Output file", required=True)
    argv = vars(parser.parse_args())
    return argv


def func(variantXls, depth, alt, freq, outf):
    # pydir = os.path.dirname(os.path.realpath(__file__))
    with open(variantXls, "r") as f, open(outf, "w") as fi:
        for line in f:
            if re.match("Chrom", line):
                fi.write(line)
                continue
            ls = line.strip().split("\t")
            if (
                ls[8] == "synonymous SNV"
                and re.search("ogenic|HOT", ":".join(ls[-6:])) is None
                or re.search("enign|eutral", ":".join(ls[-6:]))
                and 40 < float(ls[7]) < 70
                or re.search("enign|eutral", ":".join(ls[-6:]))
                and 88 < float(ls[7])
            ):
                continue
            if (
                int(float(ls[5])) < depth
                or int(float(ls[6])) < alt
                or float(ls[7]) < freq
            ):
                continue
            fi.write(line)

    # time_end = datetime.datetime.now()
    # print(
    #     f"ONE variant QC filter used: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    # )


def main():
    argv = argparse_line()
    func(
        argv["snvIndels"],
        argv["depth"],
        argv["alt"],
        argv["frequency"],
        argv["out"],
    )


if __name__ == "__main__":
    main()
