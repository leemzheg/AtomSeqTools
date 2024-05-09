#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/03/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os
import argparse, textwrap
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
import datetime, subprocess, pysam

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawTextHelpFormatter
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument("-bed", metavar="STR", help="i", required=True)
    required.add_argument("-outdir", metavar="PATH", help="o", required=True)
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "--fq-prefix",
        metavar="STR",
        help="  Specify fastq prefix to be analyzed [default: XXX]"
        "  [eg: '0101-XXX-M3-A1_1.fq.gz' prefix is '0101-XXX-M3-A1']",
        default="XXX",
    )
    argv = vars(parser.parse_args())
    return argv


def FqSummary(outd, samplename):
    stat_dic = {}
    stat_dic["Sample"] = samplename
    try:
        qc_df = pd.read_json(f"{outd}/01_CleanFq.json")
    except:
        print("Wrong: Can`t find Sample-Dir/Mid-files/01_CleanFq.json. ")
        exit()
    raw_reads_num = int(qc_df["read1_before_filtering"]["total_reads"])
    clean_reads_num = int(qc_df["read1_after_filtering"]["total_reads"])
    stat_dic["Raw(G)"] = raw_G = round(raw_reads_num * 300 / 1000000000, 2)
    if raw_G > 3:
        print(
            f"Warning: this sample size is too big ({raw_G} G), "
            f"may take a long time to calculate the depth. Please be patient and wait"
        )
    stat_dic["Raw read(PE)"] = raw_reads_num
    short_read = int(qc_df["filtering_result"]["too_short_reads"])
    total_read = int(qc_df["summary"]["before_filtering"]["total_reads"])
    stat_dic["< 50bp read%"] = round(float(short_read / total_read) * 100, 2)
    stat_dic["Clean read%"] = round(clean_reads_num / raw_reads_num * 100, 2)
    q30 = qc_df["summary"]["before_filtering"]["q30_rate"]
    stat_dic["Q30 rate%"] = round(float(q30) * 100, 2)
    return stat_dic, clean_reads_num


def BAMsummary(stat_dic, clean_reads_num, outd):
    std = subprocess.run(
        f"samtools flagstat -@ 2 {outd}/03_bwa_sort.bam|grep 'properly paired ('|awk '{{print $1}}' ",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    mapping_reads_num = int(std.stdout.strip()) / 2
    stat_dic["Mapping rate%"] = round(mapping_reads_num / clean_reads_num * 100, 2)
    std = subprocess.run(
        f"samtools flagstat -@ 2 {outd}/05_deduplicated_sort.bam|grep 'properly paired ('|awk '{{print $1}}' ",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    stat_dic["Unique read(PE)"] = unique_read = int(int(std.stdout.strip()) / 2)
    stat_dic["Unique percent%"] = round(unique_read / mapping_reads_num * 100, 2)

    return stat_dic, unique_read, mapping_reads_num


def count_match_length(read_strand, cigartuples):
    cigar_dic = {0: "M", 1: "I", 2: "D", 4: "S"}
    match_len = 0
    if read_strand == "+":
        if cigartuples[0][0] == 0 or cigartuples[0][0] == 4 and cigartuples[0][1] <= 2:
            for tu in cigartuples:
                if tu[0] == 0 or tu[0] == 2:
                    match_len += tu[1]
    else:
        if (
            cigartuples[-1][0] == 0
            or cigartuples[-1][0] == 4
            and cigartuples[-1][1] <= 2
        ):
            for tu in cigartuples:
                if tu[0] == 0 or tu[0] == 2:
                    match_len += tu[1]
    return match_len


def CountPrimerDepth(stat_dic, outd, bed, unique_read, mapping_reads_num):
    # deduplicated bam
    dic = {}
    with open(bed, "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            dic[ke] = {}
            dic[ke]["chr"] = ls[0]
            dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            dic[ke]["sta"] = int(ls[1])
            dic[ke]["end"] = int(ls[2])
            dic[ke]["primer"] = ls[3]
            dic[ke]["length"] = int(ls[4])
            dic[ke]["strand"] = ls[5]
            dic[ke]["depth"] = 0
    subprocess.run(
        f"samtools view -@ 2 -b -f 64 -F 2048 {outd}/05_deduplicated_sort.bam >{outd}/zz.tmp.sort.bam",
        shell=True,
    )
    subprocess.run(f"samtools index -@ 2 {outd}/zz.tmp.sort.bam", shell=True)
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/zz.tmp.sort.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in dic.keys():
                    if (
                        dic[k]["chr"] == chr_name
                        and dic[k]["strand"] == read_strand
                        and np.abs(position - dic[k]["pos"]) <= 2
                        and match_len > dic[k]["length"]
                    ):
                        dic[k]["depth"] += 1
                        fit_read += 1
    stat_dic["Unique on target rate%"] = round(fit_read / unique_read * 100, 2)
    depth_dic, fit = {}, 0
    for k in dic.keys():
        if dic[k]["primer"] not in depth_dic.keys():
            depth_dic[dic[k]["primer"]] = {}
            depth_dic[dic[k]["primer"]]["Unique depth(X)"] = dic[k]["depth"]
        else:
            depth_dic[dic[k]["primer"]]["Unique depth(X)"] += dic[k]["depth"]
    uniqueMeanDepth = int(fit_read / len(depth_dic))
    for k in depth_dic.keys():
        if depth_dic[k]["Unique depth(X)"] >= uniqueMeanDepth * 0.2:
            fit += 1
    UniqueUniformity = round(fit / len(depth_dic) * 100, 2)
    # bwa sort bam
    dic = {}
    with open(bed, "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            dic[ke] = {}
            dic[ke]["chr"] = ls[0]
            dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            dic[ke]["sta"] = int(ls[1])
            dic[ke]["end"] = int(ls[2])
            dic[ke]["primer"] = ls[3]
            dic[ke]["length"] = int(ls[4])
            dic[ke]["strand"] = ls[5]
            dic[ke]["depth"] = 0
    subprocess.run(
        f"samtools view -@ 2 -b -f 64 -F 2048 {outd}/03_bwa_sort.bam >{outd}/zz.tmp.sort.bam",
        shell=True,
    )
    subprocess.run(f"samtools index -@ 2 {outd}/zz.tmp.sort.bam", shell=True)
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/zz.tmp.sort.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in dic.keys():
                    if (
                        dic[k]["chr"] == chr_name
                        and dic[k]["strand"] == read_strand
                        and np.abs(position - dic[k]["pos"]) <= 2
                        and match_len > dic[k]["length"]
                    ):
                        dic[k]["depth"] += 1
                        fit_read += 1
    stat_dic["On target rate%"] = round(fit_read / mapping_reads_num * 100, 2)
    depth2_dic, fit = {}, 0
    for k in dic.keys():
        if dic[k]["primer"] not in depth2_dic.keys():
            depth2_dic[dic[k]["primer"]] = {}
            depth2_dic[dic[k]["primer"]]["Depth(X)"] = dic[k]["depth"]
        else:
            depth2_dic[dic[k]["primer"]]["Depth(X)"] += dic[k]["depth"]
    stat_dic["Unique mean depth(X)"] = uniqueMeanDepth
    stat_dic["Mean depth(X)"] = MeanDepth = int(fit_read / len(depth2_dic))
    for k in depth2_dic.keys():
        if depth2_dic[k]["Depth(X)"] >= MeanDepth * 0.2:
            fit += 1
    stat_dic["Unique uniformity%"] = UniqueUniformity
    stat_dic["Uniformity%"] = round(fit / len(depth2_dic) * 100, 2)
    unique_depth_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    depth_df = (
        pd.DataFrame(depth2_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    df = pd.merge(depth_df, unique_depth_df, on="Primer", how="outer")
    os.system(f"rm -f {outd}/zz.tmp*")
    return stat_dic, df


def func(primerbed, outd, samplename):
    # pydir = os.path.dirname(os.path.realpath(__file__))
    stat_dic, clean_reads_num = FqSummary(outd, samplename)  # method
    stat_dic, unique_read, mapping_reads_num = BAMsummary(
        stat_dic, clean_reads_num, outd
    )  # method
    stat_dic, df = CountPrimerDepth(
        stat_dic, outd, primerbed, unique_read, mapping_reads_num
    )  # method
    df.to_csv(f"{outd}/../{samplename}_primer_depth.xls", sep="\t", index=None)
    stat_df = pd.DataFrame([stat_dic])
    stat_df.to_csv(f"{outd}/../{samplename}_basic_summary.xls", sep="\t", index=None)

    time_end = datetime.datetime.now()
    print(
        f"   Read basic summary used: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    )


def main():
    argv = argparse_line()
    func(argv["bed"], argv["outdir"], argv["fq_prefix"])


if __name__ == "__main__":
    main()
