#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2023/07/25 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import sys, re, os
import argparse, textwrap
import pandas as pd


def argparse_line():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawTextHelpFormatter
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument("-i", metavar="i", help="i", required=True)
    required.add_argument("-shift", metavar="i", help="i", required=True)
    required.add_argument("-tsv", metavar="i", help="i", required=True)
    required.add_argument("-o", metavar="o", help="o", required=True)
    argv = vars(parser.parse_args())
    return argv


def func(left_vcf, shift, tsv, right_vcf):
    # pydir = os.path.dirname(os.path.realpath(__file__))
    tsv_dic = {}
    with open(tsv) as f:
        for line in f:
            if re.match("Hugo", line):
                continue
            ls = line.strip().split("\t")
            tsv_dic[ls[4]] = ""
    right_dic = {}
    with open(shift, "r") as f:
        for line in f:
            if re.match("#", line):
                continue
            ls = line.strip().split("\t")
            if (
                len(ls[-1].split(";")) < 2
                or int(ls[-1].split(";")[-1].split("=")[-1]) == 0
            ):
                continue
            else:
                if ls[4] in tsv_dic:
                    right_dic[ls[0]] = (
                        ls[1].split("-")[0]
                        + ";"
                        + ls[2]
                        + ";"
                        + ls[-3]
                        + ";"
                        + str(ls[-1].split("=")[-1])
                    )
    with open(left_vcf, "r") as f, open(right_vcf, "w") as fi:
        for line in f:
            ls = line.strip().split("\t")
            if ls[2] not in right_dic:
                fi.write("{}".format(line))
            else:
                ls[0] = right_dic[ls[2]].split(";")[0].split(":")[0]
                shift_len = int(right_dic[ls[2]].split(";")[-1])
                ref = ls[3]
                alt = ls[4]
                if len(ls[3]) > len(ls[4]):
                    ls[1] = str(int(right_dic[ls[2]].split(";")[0].split(":")[1]) - 1)
                    match = re.search(
                        r"([a-z])([A-Z]+)", right_dic[ls[2]].split(";")[-2]
                    )
                    if match:
                        ls[3] = match.group(1).upper() + match.group(2)
                        ls[4] = match.group(1).upper()
                    else:
                        if shift_len > len(ls[3]) - 1:
                            ls[3] = ref[-1] + ref[1:]
                            ls[4] = ref[-1]
                        else:
                            ls[3] = ref[shift_len:] + ref[1 : shift_len + 1]
                            ls[4] = ref[shift_len]
                    fi.write("{}\n".format("\t".join(ls)))
                elif len(ls[3]) < len(ls[4]):
                    ls[1] = right_dic[ls[2]].split(";")[0].split(":")[1]
                    match = re.search(
                        r"([a-z])([A-Z]+)", right_dic[ls[2]].split(";")[-2]
                    )
                    if match:
                        ls[3] = match.group(1).upper()
                        ls[4] = match.group(1).upper() + match.group(2).upper()
                    else:
                        if shift_len > len(ls[4]) - 1:
                            ls[3] = alt[-1]
                            ls[4] = alt[-1] + right_dic[ls[2]].split(";")[1]
                        else:
                            ls[3] = alt[shift_len]
                            ls[4] = alt[shift_len] + right_dic[ls[2]].split(";")[1]
                    fi.write("{}\n".format("\t".join(ls)))
                else:
                    fi.write("{}".format(line))


def main():
    argv = argparse_line()
    func(argv["i"], argv["shift"], argv["tsv"], argv["o"])


if __name__ == "__main__":
    main()
