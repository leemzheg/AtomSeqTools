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
import datetime, subprocess
import pysam

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawTextHelpFormatter
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-outdir", metavar="PATH", help="output directory", required=True
    )
    required.add_argument("-bed", metavar="STR", help="primer bed file", required=True)
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


def qc_from_fq(outd, samplename):
    stat_dic = {}
    stat_dic["Sample"] = samplename
    try:
        qc_df = pd.read_json(f"{outd}/01_CleanFqMeth.json")
    except:
        print("Wrong: Can`t find Sample-Dir/Mid-files/01_CleanFqMeth.json.")
        exit()
    raw_reads_num = int(qc_df["read1_before_filtering"]["total_reads"])
    clean_reads_num = int(qc_df["read1_after_filtering"]["total_reads"])
    stat_dic["Raw(G)"] = round(raw_reads_num * 300 / 1000000000, 2)
    stat_dic["Raw read(PE)"] = raw_reads_num
    short_read = int(qc_df["filtering_result"]["too_short_reads"])
    total_read = int(qc_df["summary"]["before_filtering"]["total_reads"])
    stat_dic["< 30bp read%"] = round(float(short_read / total_read) * 100, 2)
    stat_dic["Clean read%"] = round(clean_reads_num / raw_reads_num * 100, 2)
    q30 = qc_df["summary"]["before_filtering"]["q30_rate"]
    stat_dic["Q30 rate%"] = round(float(q30) * 100, 2)
    return stat_dic, clean_reads_num


def read_from_bam(outd, stat_dic, clean_reads_num):
    std = subprocess.run(
        f"samtools flagstat -@ 2 {outd}/03_bismark_sort.bam|grep 'properly paired ('|awk '{{print $1}}' ",
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


def GetMethylationDict(methy_dic, convert, tag_xm, amplicon, amplicon_dic):
    z_num, Z_num, u_num, U_num, x_num, X_num, h_num, H_num = 0, 0, 0, 0, 0, 0, 0, 0
    for ele in tag_xm:
        if ele == ".":
            continue
        elif ele == "z":
            z_num += 1
        elif ele == "Z":
            Z_num += 1
        elif ele == "u":
            u_num += 1
        elif ele == "U":
            U_num += 1
        elif ele == "x":
            x_num += 1
        elif ele == "X":
            X_num += 1
        elif ele == "h":
            h_num += 1
        elif ele == "H":
            H_num += 1
        else:
            print("tag XM had wrong")
            exit()
    methy_dic[convert]["CpG_unmethylated"] += z_num
    methy_dic[convert]["CpG_methylated"] += Z_num
    methy_dic[convert]["Unknown_unmethylated"] += u_num
    methy_dic[convert]["Unknown_methylated"] += U_num
    methy_dic[convert]["CHG_unmethylated"] += x_num
    methy_dic[convert]["CHG_methylated"] += X_num
    methy_dic[convert]["CHH_unmethylated"] += h_num
    methy_dic[convert]["CHH_methylated"] += H_num
    amplicon_dic[amplicon]["cpg_unmeth"] += z_num
    amplicon_dic[amplicon]["cpg_meth"] += Z_num
    a = u_num + x_num + h_num
    b = U_num + X_num + H_num
    amplicon_dic[amplicon]["noncpg_unmeth"] += a
    amplicon_dic[amplicon]["noncpg_meth"] += b
    cpg_converted = round(z_num / (z_num + Z_num) * 100) if (z_num + Z_num) != 0 else 0
    noncpg_converted = round(a / (a + b) * 100) if (a + b) != 0 else 0
    return methy_dic, amplicon_dic, cpg_converted, noncpg_converted


def depth_from_bam(stat_dic, outd, bed, unique_read):
    dic, methy_dic, amplicon_dic = {}, {}, {}
    convert_list = ["CT-CT", "CT-GA", "GA-GA", "GA-CT"]
    for key in convert_list:
        methy_dic[key] = {
            "CpG_unmethylated": 0,
            "CpG_methylated": 0,
            "Unknown_unmethylated": 0,
            "Unknown_methylated": 0,
            "CHG_unmethylated": 0,
            "CHG_methylated": 0,
            "CHH_unmethylated": 0,
            "CHH_methylated": 0,
        }
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
            dic[ke]["orient"] = "CT" if ls[5] == "+" else "GA"
            dic[ke]["depth"] = 0
            amplicon_dic[ls[3]] = {
                "cpg_unmeth": 0,
                "cpg_meth": 0,
                "noncpg_unmeth": 0,
                "noncpg_meth": 0,
            }
    subprocess.run(
        f"samtools view -@ 2 -b -f 64 -F 2048 {outd}/05_deduplicated_sort.bam >{outd}/zz.tmp.sort.bam",
        shell=True,
    )
    subprocess.run(f"samtools index -@ 2 {outd}/zz.tmp.sort.bam", shell=True)
    fit_read, ontarget_dic = 0, {}
    with pysam.AlignmentFile(f"{outd}/zz.tmp.sort.bam", "rb") as f:
        with pysam.AlignmentFile(f"{outd}/zz.ontarget.bam", "wb", template=f) as output:
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
                            and dic[k]["orient"] == read.get_tag("XG")
                            and np.abs(position - dic[k]["pos"]) <= 2
                            and match_len > dic[k]["length"]
                        ):
                            primername = dic[k]["primer"]
                            dic[k]["depth"] += 1
                            fit_read += 1
                            convert = read.get_tag("XR") + "-" + read.get_tag("XG")
                            methy_dic, amplicon_dic, cpg_converted, noncpg_converted = (
                                GetMethylationDict(
                                    methy_dic,
                                    convert,
                                    read.get_tag("XM"),
                                    primername,
                                    amplicon_dic,
                                )
                            )
                            ontarget_dic[fit_read] = {}
                            ontarget_dic[fit_read]["primer"] = primername
                            ontarget_dic[fit_read]["sta"] = dic[k]["sta"]
                            ontarget_dic[fit_read]["end"] = dic[k]["end"]
                            ontarget_dic[fit_read]["length"] = dic[k]["length"]
                            ontarget_dic[fit_read]["read_chr"] = chr_name
                            ontarget_dic[fit_read]["read_sta"] = chr_sta
                            ontarget_dic[fit_read]["read_end"] = chr_sta + match_len
                            ontarget_dic[fit_read]["read_cigar"] = read.cigarstring
                            ontarget_dic[fit_read]["strand"] = read_strand
                            ontarget_dic[fit_read]["primer_ref"] = dic[k]["orient"]
                            ontarget_dic[fit_read]["read_id"] = read.query_name
                            ontarget_dic[fit_read]["read_length"] = match_len
                            ontarget_dic[fit_read]["cpg_converted%"] = cpg_converted
                            ontarget_dic[fit_read][
                                "noncpg_converted%"
                            ] = noncpg_converted
                            output.write(read)
    df = pd.DataFrame(ontarget_dic).T
    df.to_csv(f"{outd}/zz.primer.ontarget.xls", sep="\t", index=None)
    stat_dic["Unique on target rate%"] = round(fit_read / unique_read * 100, 2)
    depth_dic, fit = {}, 0
    for k in dic.keys():
        if dic[k]["primer"] not in depth_dic.keys():
            depth_dic[dic[k]["primer"]] = {}
            depth_dic[dic[k]["primer"]]["Unique depth(X)"] = dic[k]["depth"]
        else:
            depth_dic[dic[k]["primer"]]["Unique depth(X)"] += dic[k]["depth"]
    uniquemeandepth = int(fit_read / len(depth_dic))
    for k in depth_dic.keys():
        if depth_dic[k]["Unique depth(X)"] >= uniquemeandepth * 0.2:
            fit += 1
        depth_dic[k]["CpG converted%"] = (
            round(
                amplicon_dic[k]["cpg_unmeth"]
                / (amplicon_dic[k]["cpg_unmeth"] + amplicon_dic[k]["cpg_meth"])
                * 100,
                2,
            )
            if (amplicon_dic[k]["cpg_unmeth"] + amplicon_dic[k]["cpg_meth"]) != 0
            else 0
        )
        depth_dic[k]["NonCpG converted%"] = (
            round(
                amplicon_dic[k]["noncpg_unmeth"]
                / (amplicon_dic[k]["noncpg_unmeth"] + amplicon_dic[k]["noncpg_meth"])
                * 100,
                2,
            )
            if (amplicon_dic[k]["noncpg_unmeth"] + amplicon_dic[k]["noncpg_meth"]) != 0
            else 0
        )
        depth_dic[k]["Conversion all%"] = (
            round(
                (amplicon_dic[k]["cpg_unmeth"] + amplicon_dic[k]["noncpg_unmeth"])
                / (
                    amplicon_dic[k]["cpg_unmeth"]
                    + amplicon_dic[k]["cpg_meth"]
                    + amplicon_dic[k]["noncpg_unmeth"]
                    + amplicon_dic[k]["noncpg_meth"]
                )
                * 100,
                2,
            )
            if (
                amplicon_dic[k]["cpg_unmeth"]
                + amplicon_dic[k]["cpg_meth"]
                + amplicon_dic[k]["noncpg_unmeth"]
                + amplicon_dic[k]["noncpg_meth"]
            )
            != 0
            else 0
        )
    UniqueUniformity = round(fit / len(depth_dic) * 100, 2)
    return stat_dic, methy_dic, depth_dic, uniquemeandepth, UniqueUniformity


def depth_from_3sort_bam(
    stat_dic, outd, bed, mapping_reads_num, uniquemeandepth, UniqueUniformity
):
    dic = {}
    with open(bed, "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            dic[ke] = {}
            dic[ke]["chr"] = ls[0]
            dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            dic[ke]["primer"] = ls[3]
            dic[ke]["length"] = int(ls[4])
            dic[ke]["strand"] = ls[5]
            dic[ke]["orient"] = "CT" if ls[5] == "+" else "GA"
            dic[ke]["depth"] = 0
    subprocess.run(
        f"samtools view -@ 2 -b -f 64 -F 2048 {outd}/03_bismark_sort.bam >{outd}/zz.tmp.sort.bam",
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
                        and dic[k]["orient"] == read.get_tag("XG")
                        and np.abs(position - dic[k]["pos"]) <= 2
                        and match_len > dic[k]["length"]
                    ):
                        dic[k]["depth"] += 1
                        fit_read += 1
    stat_dic["On target rate%"] = round(fit_read / mapping_reads_num * 100, 2)
    depth_dic, fit = {}, 0
    for k in dic.keys():
        if dic[k]["primer"] not in depth_dic.keys():
            depth_dic[dic[k]["primer"]] = {}
            depth_dic[dic[k]["primer"]]["Depth(X)"] = dic[k]["depth"]
        else:
            depth_dic[dic[k]["primer"]]["Depth(X)"] += dic[k]["depth"]
    stat_dic["Unique mean depth(X)"] = uniquemeandepth
    stat_dic["Mean depth(X)"] = MeanDepth = int(fit_read / len(depth_dic))
    for k in depth_dic.keys():
        if depth_dic[k]["Depth(X)"] >= MeanDepth * 0.2:
            fit += 1
    stat_dic["Unique uniformity%"] = UniqueUniformity
    stat_dic["Uniformity%"] = round(fit / len(depth_dic) * 100, 2)
    return stat_dic, depth_dic


def GetMethylationRate(stat_dic, v):
    methylated = {}
    # all stand
    cpg = (
        v["CT-CT"]["CpG_methylated"]
        + v["GA-GA"]["CpG_methylated"]
        + v["CT-GA"]["CpG_methylated"]
        + v["GA-CT"]["CpG_methylated"]
    )
    chg = (
        v["CT-CT"]["CHG_methylated"]
        + v["GA-GA"]["CHG_methylated"]
        + v["CT-GA"]["CHG_methylated"]
        + v["GA-CT"]["CHG_methylated"]
    )
    chh = (
        v["CT-CT"]["CHH_methylated"]
        + v["GA-GA"]["CHH_methylated"]
        + v["CT-GA"]["CHH_methylated"]
        + v["GA-CT"]["CHH_methylated"]
    )
    methylated["CpG_methylated%"] = (
        cpg
        / (
            cpg
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
        )
        * 100
    ) if (
            cpg
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
        ) != 0 else 0
    methylated["nonCpG_methylated%"] = (
        (chg + chh)
        / (
            chg
            + chh
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            chg
            + chh
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        ) != 0 else 0
    methylated["methylated_all%"] = (
        (cpg + chg + chh)
        / (
            cpg
            + chg
            + chh
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            cpg
            + chg
            + chh
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        ) != 0 else 0
    # CT GA strand
    cpg_ct = v["CT-CT"]["CpG_methylated"] + v["GA-CT"]["CpG_methylated"]
    cpg_ga = v["GA-GA"]["CpG_methylated"] + v["CT-GA"]["CpG_methylated"]
    chg_ct = v["CT-CT"]["CHG_methylated"] + v["GA-CT"]["CHG_methylated"]
    chg_ga = v["GA-GA"]["CHG_methylated"] + v["CT-GA"]["CHG_methylated"]
    chh_ct = v["CT-CT"]["CHH_methylated"] + v["GA-CT"]["CHH_methylated"]
    chh_ga = v["GA-GA"]["CHH_methylated"] + v["CT-GA"]["CHH_methylated"]
    methylated["CpG_methylated_CT%"] = (
        cpg_ct
        / (cpg_ct + v["CT-CT"]["CpG_unmethylated"] + v["GA-CT"]["CpG_unmethylated"])
        * 100
    ) if (cpg_ct + v["CT-CT"]["CpG_unmethylated"] + v["GA-CT"]["CpG_unmethylated"]) != 0 else 0
    methylated["CpG_methylated_GA%"] = (
        cpg_ga
        / (cpg_ga + v["CT-GA"]["CpG_unmethylated"] + v["GA-GA"]["CpG_unmethylated"])
        * 100
    ) if (cpg_ga + v["CT-GA"]["CpG_unmethylated"] + v["GA-GA"]["CpG_unmethylated"]) != 0 else 0
    methylated["nonCpG_methylated_CT%"] = (
        (chg_ct + chh_ct)
        / (
            chg_ct
            + chh_ct
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            chg_ct
            + chh_ct
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        ) != 0 else 0
    methylated["nonCpG_methylated_GA%"] = (
        (chg_ga + chh_ga)
        / (
            chg_ga
            + chh_ga
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            chg_ga
            + chh_ga
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
        ) != 0 else 0
    methylated["Methylation_CT%"] = (
        (cpg_ct + chg_ct + chh_ct)
        / (
            cpg_ct
            + chg_ct
            + chh_ct
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            cpg_ct
            + chg_ct
            + chh_ct
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        ) != 0 else 0
    methylated["Methylation_GA%"] = (
        (cpg_ga + chg_ga + chh_ga)
        / (
            cpg_ga
            + chg_ga
            + chh_ga
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
        )
        * 100
    ) if (
            cpg_ga
            + chg_ga
            + chh_ga
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
        ) != 0 else 0
    stat_dic["CT strand CpG converted%"] = round(
        100 - methylated["CpG_methylated_CT%"], 2
    )
    stat_dic["GA strand CpG converted%"] = round(
        100 - methylated["CpG_methylated_GA%"], 2
    )
    stat_dic["Both strand CpG converted%"] = round(
        100 - methylated["CpG_methylated%"], 2
    )
    stat_dic["CT strand nonCpG converted%"] = round(
        100 - methylated["nonCpG_methylated_CT%"], 2
    )
    stat_dic["GA strand nonCpG converted%"] = round(
        100 - methylated["nonCpG_methylated_GA%"], 2
    )
    stat_dic["Both strand nonCpG converted%"] = round(
        100 - methylated["nonCpG_methylated%"], 2
    )
    stat_dic["CT strand all converted%"] = round(100 - methylated["Methylation_CT%"], 2)
    stat_dic["GA strand all converted%"] = round(100 - methylated["Methylation_GA%"], 2)
    stat_dic["Both strand all converted%"] = round(
        100 - methylated["methylated_all%"], 2
    )
    return stat_dic


def func(outd, samplename, bed):
    # pydir = os.path.dirname(os.path.realpath(__file__))

    stat_dic, clean_reads_num = qc_from_fq(outd, samplename)  # method
    stat_dic, unique_read, mapping_reads_num = read_from_bam(
        outd, stat_dic, clean_reads_num
    )  # method
    stat_dic, methy_dic, depth_dic, uniquemeandepth, UniqueUniformity = depth_from_bam(
        stat_dic, outd, bed, unique_read
    )  # method

    depth_5dedup_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    stat_dic, depth_dic = depth_from_3sort_bam(
        stat_dic,
        outd,
        bed,
        mapping_reads_num,
        uniquemeandepth,
        UniqueUniformity,
    )  # method
    depth_3sort_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    depth_df = pd.merge(depth_3sort_df, depth_5dedup_df, on="Primer", how="outer")
    depth_df["Unique depth(X)"] = depth_df["Unique depth(X)"].astype(int)
    depth_df = depth_df.sort_values(by="Primer", ascending=True)
    depth_df.to_csv(f"{outd}/../{samplename}_depth.xls", sep="\t", index=None)

    stat_dic = GetMethylationRate(stat_dic, methy_dic)  # method
    stat_df = pd.DataFrame([stat_dic])
    stat_df.to_csv(f"{outd}/../{samplename}_basic_summary.xls", sep="\t", index=None)
    os.system(f"rm -f {outd}/zz.* ")

    time_end = datetime.datetime.now()
    print(
        f"   ONE methylation summary used: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    )


def main():
    argv = argparse_line()
    func(
        argv["outdir"],
        argv["fq_prefix"],
        argv["bed"],
    )


if __name__ == "__main__":
    main()
