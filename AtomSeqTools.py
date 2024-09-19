#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/08/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os, yaml
import argparse, textwrap
import datetime, glob
from pathlib import Path
import subprocess, gzip
import pandas as pd
import numpy as np
import pysam

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="Run AtomSeq Target Pipeline Toolkit",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="The bed file of primers should contain at least first six columns(delimited by table):\n"
        "1. Chormosome ID  (consistent with reference genome);\n"
        "2. Start of primer (0-based);\n"
        "3. End of primer (1-based);\n"
        "4. Primer ID;\n"
        "5. Primer length;\n"
        "6. Targeted strand (+/-).",
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-image", metavar="IMAGE", help="    Singularity image [*.sif]", required=True
    )
    required.add_argument(
        "-fq-dir",
        metavar="PATH",
        help="    The path to directory of raw fastq",
        required=True,
    )
    required.add_argument(
        "-fq-prefix",
        metavar="STR",
        help="    Specify fastq prefix to be analyzed (can be many, delimited by space.\n"
        "    If baselineCreator mode, the more healthy samples, the better)\n"
        "    [eg: '0101-XXX-M3-A1_1.fq.gz' prefix is '0101-XXX-M3-A1']",
        required=True,
        nargs="+",
    )
    required.add_argument(
        "-outdir",
        metavar="PATH",
        help="    The path to directory of all output files",
        required=True,
    )
    required.add_argument(
        "-config",
        metavar="STR",
        help="    Config file. Some parameters and paths that need to be set",
        required=True,
    )
    required.add_argument(
        "-bed",
        metavar="STR",
        help="    Input bed file of primers",
        required=True,
    )
    variant = parser.add_argument_group("Variant options")
    variant.add_argument(
        "-variant",
        help="    Enable pipeline to call variant",
        action="store_true",
    )
    variant.add_argument(
        "-sample-type",
        metavar="STR",
        help="    If call variant, it`s necessary to select sample type. Possible values:{cfDNA, Tissue}",
        choices=["cfDNA", "Tissue"],
        required=False,
    )
    variant.add_argument(
        "-cnv",
        action="store_true",
        help="    Enable pipeline to call cnv (cfDNA can`t call cnv)",
    )
    variant.add_argument(
        "-cnv-baseline",
        metavar="STR",
        help="    If call cnv, it`s best to give a cnv baseline, but can be None. [default: None]",
        required=False,
        default=None,
    )
    variant.add_argument(
        "-msi",
        action="store_true",
        help="    Enable pipeline to call msi (cfDNA can`t call msi)",
    )
    variant.add_argument(
        "-msi-baseline",
        metavar="STR",
        help="    If call msi, it must give a msi baseline",
        required=False,
    )
    meth = parser.add_argument_group("Methylation options")
    meth.add_argument(
        "-meth",
        action="store_true",
        help="    Enable pipeline to analyse methylation data",
    )
    fusion = parser.add_argument_group("Fusion options")
    fusion.add_argument(
        "-fusion",
        action="store_true",
        help="    Enable pipeline to analyse RNA fusion data(When running this mode, 40G of \n"
        "    computer memory is required, otherwise it will crash and interrupt)",
    )
    make_baseline = parser.add_argument_group("BaselineCreator options")
    make_baseline.add_argument(
        "-baselineCreator",
        action="store_true",
        help="    Enable pipeline to make a CNV/MSI baseline file",
    )
    make_baseline.add_argument(
        "-baseline-type",
        metavar="STR",
        help="    Select output baseline type. Possible values:{CNV, MSI}",
        choices=["CNV", "MSI"],
        required=False,
    )
    make_baseline.add_argument(
        "-baseline-name",
        metavar="STR",
        help="    The name of the CNV/MSI baseline [default: XXXBaseline_XXXpanel_AtomSeq.XXXgene]",
        default="Baseline_XXXpanel_AtomSeq.XXXgene",
    )
    optional = parser.add_argument_group("Other options")
    optional.add_argument(
        "-threads",
        metavar="INT",
        help="    Number of threads to use [default:8]",
        type=int,
        default=8,
    )
    optional.add_argument(
        "-data-to-analyse",
        metavar="FLOAT",
        help="    Specify how many data size(G) to be analysed \n"
        "    [default:0] [0 means analyse all fastq data]",
        type=float,
        default=0,
    )
    optional.add_argument(
        "-supporting-reads",
        metavar="INT",
        help="    Only output consensus reads/pairs that merged by >= <supporting-reads> reads/pairs.\n"
        "    The valud should be 1~10, and the default value is 1. [default: 1]",
        type=int,
        default=1,
    )
    optional.add_argument(
        "-auto-remove",
        action="store_true",
        help="    Remove all temporary files (will remove Mid-files)",
    )
    argv = vars(parser.parse_args())
    return argv


def GetConfigFile(config, mount_paras, image, pydir):
    config_dic = {}
    with open(config, "r") as f:
        for line in f:
            if re.match("#", line):
                continue
            if re.search("=", line):
                ls = line.strip().split("=")
                config_dic[ls[0]] = ls[1]
                if ls[0] == "Hg38_Fasta_Path":
                    # tmp_path = "/".join(ls[1].split("/")[:-1])
                    tmp_path = os.path.dirname(ls[1])
                    mount_paras += f" -B {tmp_path}"
                    bwa_dir = Path(f"{tmp_path}/Bwa_index")
                    bismark_dir = Path(f"{tmp_path}/Bisulfite_Genome")
                    if not bwa_dir.exists() or not bismark_dir.exists():
                        if not bwa_dir.exists():
                            bwa_dir.mkdir(parents=True)
                        cmd = (
                            f"python3 {pydir}/bin/make_index.py "
                            f"-image {image} -fasta {ls[1]}"
                        )
                        print(
                            f"   {cmd}\n\033[33mNotice\033[0m: Detected that hg38 "
                            "alignment index doesn't exist. Building for you, "
                            "taking approximately 100 minutes "
                        )
                        os.system(cmd)
                else:
                    mount_paras += f" -B {ls[1]}"
    return config_dic, mount_paras


def FastqDispose(fqdir, sample_name):
    sample_info = {}
    for fq in glob.glob(f"{fqdir}/*_*1.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        suffix = "_" + "_".join(fq.split("/")[-1].split("_")[1:])
        if sample not in sample_info:
            sample_info[sample] = {}
            sample_info[sample]["directory"] = fqdir
            sample_info[sample]["prefix"] = sample
            sample_info[sample]["R1_suffix"] = suffix
        else:
            print(
                f"The {sample} sample name in the rawdata directory is duplicated, please check"
            )
            continue
    for fq in glob.glob(f"{fqdir}/*_*2.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        suffix = "_" + "_".join(fq.split("/")[-1].split("_")[1:])
        if sample not in sample_info:
            print(
                f"In the rawdata directory, {sample} only has fq2 but not fq1. Please check"
            )
            exit()
        else:
            sample_info[sample]["R2_suffix"] = suffix
    non_saved = []
    if sample_name:
        for name in sample_name:
            if name not in list(sample_info.keys()):
                non_saved.append(name)
                continue
        if len(non_saved) > 0:
            print(
                f"Traceback: fastq directory don`t have {'、'.join(non_saved)}, please check"
            )
            exit()
        for sample in list(sample_info.keys()):
            if sample not in sample_name:
                del sample_info[sample]
                continue
    # 检出字典是否为空
    if not sample_info:
        print(
            f"{fqdir} fastq format has wrong or {'、'.join(sample_name)} has wrong, please check"
        )
        exit()
    return sample_info


def CheckModeConflict(
    variantCall,
    cnvCall,
    msiCall,
    methCall,
    fusionCall,
    sampleType,
    msiBaseline,
    baselineCall,
    baselineType,
):
    if variantCall or cnvCall or msiCall:
        if variantCall:
            if not sampleType:
                print(
                    "Traceback: when variant calling, -sample-type is required. "
                    "Select sample type, input 'cfDNA' or 'Tissue'. "
                    r"Possible values: {cfDNA, Tissue}"
                )
                exit()
        if cnvCall:
            if sampleType == "cfDNA":
                print("Warning: Not recommended to call CNV for cfDNA")
        if msiCall:
            if not msiBaseline:
                print(
                    "Traceback: when call msi, -msi-baseline is required. "
                    "Please give a msi baseline"
                )
                exit()
            if sampleType == "cfDNA":
                print("Warning: Not recommended to call MSI for cfDNA")
        if methCall:
            print(
                "Traceback: You have already opened Variant/cnv/msi mode, you can`t open "
                "Methylation mode (-meth) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if fusionCall:
            print(
                "Traceback: You have already opened Variant/cnv/msi mode, you can`t open "
                "Fusion mode (-fusion) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineCall:
            print(
                "Traceback: You have already opened Variant/cnv/msi mode, you can`t open "
                "baselineCreator mode (-baselineCreator) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineType:
            print(
                "Traceback: You have already opened Variant/cnv/msi mode, you can`t open "
                "baselineCreator mode (-baseline-type) at the same time(Only one mode can be activated), please check"
            )
            exit()
    if methCall:
        if fusionCall:
            print(
                "Traceback: You have already opened Methylation mode (-meth), you can`t open "
                "Fusion mode (-fusion) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if variantCall or cnvCall or msiCall:
            print(
                "Traceback: You have already opened Methylation mode (-meth), you can`t open "
                "Variant/cnv/msi mode at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineCall:
            print(
                "Traceback: You have already opened Methylation mode (-meth), you can`t open "
                "baselineCreator mode (-baselineCreator) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineType:
            print(
                "Traceback: You have already opened Methylation mode, you can`t open "
                "baselineCreator mode (-baseline-type) at the same time(Only one mode can be activated), please check"
            )
            exit()
    if fusionCall:
        if methCall:
            print(
                "Traceback: You have already opened Fusion mode (-fusion), you can`t open "
                "Methylation mode (-meth) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if variantCall or cnvCall or msiCall:
            print(
                "Traceback: You have already opened Fusion mode (-fusion), you can`t open "
                "Variant/cnv/msi mode at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineCall:
            print(
                "Traceback: You have already opened Fusion mode (-fusion), you can`t open "
                "baselineCreator mode (-baselineCreator) at the same time(Only one mode can be activated), please check"
            )
            exit()
        if baselineType:
            print(
                "Traceback: You have already opened Fusion mode (-fusion), you can`t open "
                "baselineCreator mode (-baseline-type) at the same time(Only one mode can be activated), please check"
            )
            exit()
    if baselineCall or baselineType:
        if baselineCall and not baselineType:
            print(
                "Traceback: You have already opened baselineCreator mode (-baselineCreator), you have to add -baseline-type "
                "Possible values:{CNV, MSI} , please check"
            )
            exit()
        if baselineType and not baselineCall:
            print(
                "Warning: You add -baseline-type but not enable pipeline to make a CNV/MSI baseline file, please "
                "add -baselineCreator"
            )
            exit()


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


def summary_variant(outd, dic, sample, mount_paras):
    # fq aspect
    stat_dic = {}
    try:
        qc_df = pd.read_json(f"{outd}/Mid-files/01_CleanFq.json")
    except:
        print("Wrong: Can`t find Sample-Dir/Mid-files/01_CleanFq.json. ")
        exit()
    raw_reads_num = int(qc_df["read1_before_filtering"]["total_reads"])
    clean_reads_num = int(qc_df["read1_after_filtering"]["total_reads"])
    raw_G = round(raw_reads_num * 300 / 1000000000, 2)
    if raw_G > 3:
        print(
            f"Warning: this sample size is too big ({raw_G} G), "
            f"may take a long time to calculate the depth. Please be patient and wait"
        )
    short_read = int(qc_df["filtering_result"]["too_short_reads"])
    total_read = int(qc_df["summary"]["before_filtering"]["total_reads"])
    short_rate = round(float(short_read / total_read) * 100, 2)
    clean_rate = round(clean_reads_num / raw_reads_num * 100, 2)
    q30 = qc_df["summary"]["before_filtering"]["q30_rate"]
    # bwa aspect
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools flagstat -@ 2 {outd}/Mid-files/03_bwa_sort.bam >{outd}/Mid-files/zz.stat.xls1"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools flagstat -@ 2 {outd}/Mid-files/05_dedup_sort.bam >{outd}/Mid-files/zz.stat.xls2"
    )
    std = subprocess.run(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"""awk '$0~"properly paired"{{print $1}}' {outd}/Mid-files/zz.stat.xls1""",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    mapping_reads_num = int(std.stdout.strip()) / 2
    map_rate = round(mapping_reads_num / clean_reads_num * 100, 2)
    std = subprocess.run(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"""awk '$0~"properly paired"{{print $1}}' {outd}/Mid-files/zz.stat.xls2""",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    unique_read = int(int(std.stdout.strip()) / 2)
    unique_percend = round(unique_read / mapping_reads_num * 100, 2)
    # primer depth dedup
    bed_dic = {}
    with open(dic["bed"], "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            bed_dic[ke] = {}
            bed_dic[ke]["chr"] = ls[0]
            bed_dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            bed_dic[ke]["sta"] = int(ls[1])
            bed_dic[ke]["end"] = int(ls[2])
            bed_dic[ke]["primer"] = ls[3]
            bed_dic[ke]["length"] = int(ls[4])
            bed_dic[ke]["strand"] = ls[5]
            bed_dic[ke]["depth"] = 0
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -@ 2 -b -f 64 -F 2048 "
        f"{outd}/Mid-files/05_dedup_sort.bam >{outd}/Mid-files/zz.stat.bam"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/zz.stat.bam"
    )
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/Mid-files/zz.stat.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in bed_dic.keys():
                    if (
                        bed_dic[k]["chr"] == chr_name
                        and bed_dic[k]["strand"] == read_strand
                        and np.abs(position - bed_dic[k]["pos"]) <= 2
                        and match_len > bed_dic[k]["length"]
                    ):
                        bed_dic[k]["depth"] += 1
                        fit_read += 1
    uni_on_target_rate = round(fit_read / unique_read * 100, 2)
    depth_dic, fit = {}, 0
    for k in bed_dic.keys():
        if bed_dic[k]["primer"] not in depth_dic.keys():
            depth_dic[bed_dic[k]["primer"]] = {}
            depth_dic[bed_dic[k]["primer"]]["Unique depth(X)"] = bed_dic[k]["depth"]
        else:
            depth_dic[bed_dic[k]["primer"]]["Unique depth(X)"] += bed_dic[k]["depth"]
    uni_mean_depth = int(fit_read / len(depth_dic))
    for k in depth_dic.keys():
        if depth_dic[k]["Unique depth(X)"] >= uni_mean_depth * 0.2:
            fit += 1
    Uni_uniformity = round(fit / len(depth_dic) * 100, 2)
    # primer depth non dedup
    bed_dic = {}
    with open(dic["bed"], "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            bed_dic[ke] = {}
            bed_dic[ke]["chr"] = ls[0]
            bed_dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            bed_dic[ke]["sta"] = int(ls[1])
            bed_dic[ke]["end"] = int(ls[2])
            bed_dic[ke]["primer"] = ls[3]
            bed_dic[ke]["length"] = int(ls[4])
            bed_dic[ke]["strand"] = ls[5]
            bed_dic[ke]["depth"] = 0
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -@ 2 -b -f 64 -F 2048 "
        f"{outd}/Mid-files/03_bwa_sort.bam >{outd}/Mid-files/zz.stat.bam"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/zz.stat.bam"
    )
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/Mid-files/zz.stat.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in bed_dic.keys():
                    if (
                        bed_dic[k]["chr"] == chr_name
                        and bed_dic[k]["strand"] == read_strand
                        and np.abs(position - bed_dic[k]["pos"]) <= 2
                        and match_len > bed_dic[k]["length"]
                    ):
                        bed_dic[k]["depth"] += 1
                        fit_read += 1
    on_target_rate = round(fit_read / mapping_reads_num * 100, 2)
    depth2_dic, fit = {}, 0
    for k in bed_dic.keys():
        if bed_dic[k]["primer"] not in depth2_dic.keys():
            depth2_dic[bed_dic[k]["primer"]] = {}
            depth2_dic[bed_dic[k]["primer"]]["Depth(X)"] = bed_dic[k]["depth"]
        else:
            depth2_dic[bed_dic[k]["primer"]]["Depth(X)"] += bed_dic[k]["depth"]
    mean_depth = int(fit_read / len(depth2_dic))
    for k in depth2_dic.keys():
        if depth2_dic[k]["Depth(X)"] >= mean_depth * 0.2:
            fit += 1
    uniformity = round(fit / len(depth2_dic) * 100, 2)
    unique_depth_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    depth_df = (
        pd.DataFrame(depth2_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    df = pd.merge(depth_df, unique_depth_df, on="Primer", how="outer")
    df.to_csv(f"{outd}/{sample}_primer_depth.xls", sep="\t", index=None)
    # final aspect
    stat_dic["Sample"] = sample
    stat_dic["Raw(G)"] = raw_G
    stat_dic["Raw read(PE)"] = raw_reads_num
    stat_dic["<50 bp read%"] = short_rate
    stat_dic["Clean read%"] = clean_rate
    stat_dic["Q30 rate%"] = round(float(q30) * 100, 2)
    stat_dic["Mapping rate%"] = map_rate
    stat_dic["Unique read(PE)"] = unique_read
    stat_dic["Unique percent%"] = unique_percend
    stat_dic["Unique on target rate%"] = uni_on_target_rate
    stat_dic["On target rate%"] = on_target_rate
    stat_dic["Unique mean depth(X)"] = uni_mean_depth
    stat_dic["Mean depth(X)"] = mean_depth
    stat_dic["Unique uniformity%"] = Uni_uniformity
    stat_dic["Uniformity%"] = uniformity
    stat_df = pd.DataFrame([stat_dic])
    stat_df.to_csv(f"{outd}/{sample}_basic_summary.xls", sep="\t", index=None)


def variant_1to5(outd, dic, sample, mount_paras):
    # fastp
    fq1 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R1_suffix"]}'
    fq2 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R2_suffix"]}'
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"fastp -i {fq1} -I {fq2} -o {outd}/Mid-files/01_CleanFq_1.fq.gz "
        f"-O {outd}/Mid-files/01_CleanFq_2.fq.gz "
        f"--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        f"--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--correction --fix_mgi_id --cut_tail -l 50 -w {dic['fastp']['threads']} "
        f"--overlap_len_require 10 "
        f"--reads_to_process {dic['fastp']['read_to_process']} -j {outd}/Mid-files/01_CleanFq.json "
        f"-h {outd}/Mid-files/01_CleanFq.html 2> {outd}/logs/01_fastp.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # add umi
    counta = 0
    with gzip.open(
        f"{outd}/Mid-files/02_UmiFq_1.fq.gz", "wt", compresslevel=4
    ) as fq1_output:
        with gzip.open(
            f"{outd}/Mid-files/02_UmiFq_2.fq.gz", "wt", compresslevel=4
        ) as fq2_output:
            with gzip.open(f"{outd}/Mid-files/01_CleanFq_1.fq.gz", "rt") as handle1:
                with gzip.open(f"{outd}/Mid-files/01_CleanFq_2.fq.gz", "rt") as handle2:
                    for tu in zip(handle1, handle2):
                        counta += 1
                        if counta % 4 == 1:
                            title1 = tu[0].split()[0]
                            title2 = tu[1].split()[0]
                        if counta % 4 == 2:
                            tmp_umi = ":UMI_" + tu[1][:20]
                            fq1_output.write(title1 + tmp_umi + "\n")
                            fq1_output.write(tu[0])
                            fq2_output.write(title2 + tmp_umi + "\n")
                            fq2_output.write(tu[1])
                        if counta % 4 == 3:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
                        if counta % 4 == 0:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
    # bwa
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bwa mem -K 10000000 -t {dic['bwa']['threads']} "
        f"-R {dic['bwa']['read_group']} {dic['bwa']['bwa_index']} "
        f"{outd}/Mid-files/02_UmiFq_1.fq.gz {outd}/Mid-files/02_UmiFq_2.fq.gz "
        f"-o {outd}/Mid-files/03_temp.sam 2> {outd}/logs/03_bwa.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # samtools sort
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -bh -f 3 -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/03_temp.sam -o {outd}/Mid-files/03_bwa_filter.bam "
    )
    # print(cmd)
    os.system(cmd)
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools sort -@ {dic['samtools']['threads']} {outd}/Mid-files/03_bwa_filter.bam "
        f"-o {outd}/Mid-files/03_bwa_sort.bam 2> {outd}/logs/03_bwa_sort.stderr "
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/03_bwa_sort.bam "
    )
    # gatk recall
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"gatk BaseRecalibrator -R {dic['bwa']['fasta']} "
        f"--known-sites {dic['gatk']['known_snps']} "
        f"--known-sites {dic['gatk']['known_indels']} "
        f"-I {outd}/Mid-files/03_bwa_sort.bam -O {outd}/Mid-files/04_recal_table.tsv "
        f"1> {outd}/logs/04_gatk_recall.stdout 2> {outd}/logs/04_gatk_recall.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # gatk bqsr
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"gatk ApplyBQSR -R {dic['bwa']['fasta']} "
        f"-bqsr {outd}/Mid-files/04_recal_table.tsv "
        f"-I {outd}/Mid-files/03_bwa_sort.bam -O {outd}/Mid-files/04_bqsr.bam "
        f"2> {outd}/logs/04_gatk_bqsr.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # gencore
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"gencore -r {dic['bwa']['fasta']} "
        f"-u UMI -s {dic['gencore']['support_reads']} "
        f"-i {outd}/Mid-files/04_bqsr.bam -o {outd}/Mid-files/05_dedup_tmp.bam "
        f"-h {outd}/Mid-files/05_dedup.html -j {outd}/Mid-files/05_dedup.json "
        f"2> {outd}/logs/05_gencore.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # sort
    if dic["baseline-creator"]:
        # samtools sort
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"samtools sort -@ {dic['samtools']['threads']} {outd}/Mid-files/05_dedup_tmp.bam "
            f"-o {outd}/Mid-files/{sample}_dedup_sort.bam 2> {outd}/logs/05_samtools_sort.stderr "
        )
        # print(cmd)
        os.system(cmd)
        # index
        os.system(
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"samtools index {outd}/Mid-files/{sample}_dedup_sort.bam "
        )
    else:
        # samtools sort
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"samtools sort -@ {dic['samtools']['threads']} {outd}/Mid-files/05_dedup_tmp.bam "
            f"-o {outd}/Mid-files/05_dedup_sort.bam 2> {outd}/logs/05_samtools_sort.stderr "
        )
        # print(cmd)
        os.system(cmd)
        # index
        os.system(
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"samtools index {outd}/Mid-files/05_dedup_sort.bam "
        )
    os.system(f"rm -f {outd}/Mid-files/03_temp.sam {outd}/Mid-files/05_dedup_tmp.bam")


def primer_convert(outd, dic, mount_paras):
    with open(dic["bed"]) as f, open(f"{outd}/zz.region.xls", "w") as fi:
        for line in f:
            ls = line.strip().split("\t")
            if ls[5] == "+":
                fi.write("{}\t{}\t{}\n".format(ls[0], ls[2], str(int(ls[2]) + 80)))
            elif ls[5] == "-":
                fi.write("{}\t{}\t{}\n".format(ls[0], str(int(ls[1]) - 80), ls[1]))
            else:
                print("Wrong: BED format has wrong, please check")
                exit()
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bedtools sort -i {outd}/zz.region.xls >{outd}/zz.region.xls2"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bedtools merge -i {outd}/zz.region.xls2 > {outd}/zz.region.bed"
    )


def variant_cnv(outd, dic, sample, mount_paras):
    if dic["cnv-baseline"]:
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"cnvkit.py batch {outd}/Mid-files/05_dedup_sort.bam "
            f"-m amplicon -r {dic['cnv-baseline']} "
            f"-d {outd}/CNVcalling 2> {outd}/logs/06_cnv.stderr"
        )
        # print(cmd)
        os.system(cmd)
    else:
        os.system(
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"cnvkit.py target {outd}/../zz.region.bed --annotate "
            f"{dic['docs']}/refFlat.txt -o {outd}/zz.tmp.bait "
            f"2> {outd}/logs/06_cnv.stderr"
        )
        os.system(
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"python3 {dic['scripts']}/mark_bed_cnvfusion.py --input "
            f"{outd}/zz.tmp.bait --cnv-gene {dic['docs']}/general_CNVgene.txt "
            f"--output {outd}/zz.tmp.bait.bed 2>> {outd}/logs/06_cnv.stderr"
        )
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"cnvkit.py batch {outd}/Mid-files/05_dedup_sort.bam "
            f"-m amplicon -f {dic['bwa']['fasta']} -t {outd}/zz.tmp.bait.bed "
            f"-n --target-avg-size 50 -d {outd}/CNVcalling 2>> {outd}/logs/06_cnv.stderr"
        )
        # print(cmd)
        os.system(cmd)
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"python3 {dic['scripts']}/calculate_cnv.py -cnr {outd}/CNVcalling/*cnr "
        f"-cnv {outd}/zz.cnv.xls "
        f"-png {outd}/{sample}_CNV_chart.png 2>> {outd}/logs/06_cnv.stderr"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"awk 'NR==1;NR>1&&$2>3.5&&$3>3.5||NR>1&&$2<0.6&&$3<0.6' {outd}/zz.cnv.xls "
        f"> {outd}/{sample}_CNV_summary.xls"
    )


def variant_msi(outd, dic, sample, mount_paras):
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"msisensor-pro pro -c 15 -d {dic['msi-baseline']} "
        f"-t {outd}/Mid-files/05_dedup_sort.bam -o {outd}/{sample}_MSI_summary.xls "
        f"1> {outd}/logs/06_msi.stdout 2> {outd}/logs/06_msi.stderr"
    )
    # print(cmd)
    os.system(cmd)
    os.system(f"rm -f {outd}/{sample}_MSI_summary.xls_*")


def variant_6to11(outd, dic, sample, mount_paras):
    # primer trim
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['scripts']}/AtomSeqTrimmer_v1.0.2.pl -primer {dic['bed']} "
        f"-thread {dic['samtools']['threads']} -type pe -split 1 -max_mis 5 "
        f"-in {outd}/Mid-files/05_dedup_sort.bam -out "
        f"{outd}/Mid-files/06_trimPrimer.bam 2> {outd}/logs/06_trimmer.stderr "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools merge -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/06_trimPrimer_temp.bam {outd}/Mid-files/06_trimPrimer_plus.bam "
        f"{outd}/Mid-files/06_trimPrimer_minus.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools sort -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/06_trimPrimer_temp.bam > {outd}/Mid-files/06_trimPrimer.bam "
        f"2> {outd}/logs/06_trimPrimer.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # varscan call vcf
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -@ {dic['samtools']['threads']} -b -q 10  "
        f"{outd}/Mid-files/06_trimPrimer.bam >{outd}/Mid-files/06_trimPrimer_q10.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools mpileup -B --min-BQ 0 --max-depth 30000 "
        f"-f {dic['bwa']['fasta']} {outd}/Mid-files/06_trimPrimer_q10.bam "
        f"-o {outd}/Mid-files/06_trimPrimer_mpileup "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"varscan mpileup2cns {outd}/Mid-files/06_trimPrimer_mpileup "
        f"--min-coverage {dic['varscan2']['min_coverage']} --min-reads2 1 "
        f"--min-avg-qual {dic['varscan2']['min_ave_quality']} --min-var-freq "
        f"0.0001 --p-value 1 --strand-filter 1 --output-vcf 1 "
        f"--variants > {outd}/Mid-files/07_variant.vcf 2> {outd}/logs/07_vcf.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # shift vcf
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bedtools intersect -a {outd}/Mid-files/07_variant.vcf "
        f"-b {outd}/../zz.region.bed -wa > {outd}/Mid-files/07_variant_tmp.vcf "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"""awk -v OFS="\\t" '{{$3="V"NR;$6="V"NR;print $0}}' """
        f"{outd}/Mid-files/07_variant_tmp.vcf > {outd}/Mid-files/07_variant_left.vcf "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"vep --species homo_sapiens --assembly GRCh38 "
        f"--hgvs --format vcf --offline --shift_3prime 1 "
        f"--force_overwrite --shift_length --shift_genomic "
        f"--dir {dic['vep']['vep_library']} --fasta {dic['bwa']['fasta']} "
        f"--input_file {outd}/Mid-files/07_variant_left.vcf "
        f"--output_file {outd}/Mid-files/07_variant_shift.vcf "
        f"1> {outd}/logs/07_vep.stdout 2> {outd}/logs/07_vep.stderr "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"python3 {dic['scripts']}/turn_vcf_to_right.py -tsv {dic['docs']}/cancerGeneList.tsv "
        f"-shift {outd}/Mid-files/07_variant_shift.vcf -i {outd}/Mid-files/07_variant_left.vcf "
        f"-o {outd}/Mid-files/07_variant_right.vcf "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"java -Xmx32g -jar {dic['snpEff']['path']} -v GRCh38.p14 -noStats "
        f"{outd}/Mid-files/07_variant_right.vcf > {outd}/Mid-files/07_variant.snpeff.vcf "
        f"2> {outd}/logs/07_snpeff.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # avinput fix
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/convert2annovar.pl -format vcf4  "
        f"{outd}/Mid-files/07_variant_left.vcf > {outd}/Mid-files/08_avinput_left.xls "
        f"2> {outd}/logs/08_avinput.stderr "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/convert2annovar.pl -format vcf4  "
        f"{outd}/Mid-files/07_variant_right.vcf > {outd}/Mid-files/08_avinput_right.xls "
        f"2>> {outd}/logs/08_avinput.stderr "
    )
    # print(cmd)
    os.system(cmd)
    dic_fix = {}
    with open(f"{outd}/Mid-files/07_variant_right.vcf") as f:
        for line in f:
            ls = line.strip().split("\t")
            dic_fix[ls[5]] = (
                ls[-1].split(":")[3]
                + ":"
                + ls[-1].split(":")[5]
                + ":"
                + ls[-1].split(":")[6].split("%")[0]
            )
    with open(f"{outd}/Mid-files/08_avinput_right.xls") as f, open(
        f"{outd}/Mid-files/08_avinput_fix.xls", "w"
    ) as fi:
        for line in f:
            ls = line.strip().split("\t")
            if ls[6] in dic_fix:
                fi.write(
                    "{}\t{}\t{}\n".format(
                        "\t".join(ls[:5]), "\t".join(dic_fix[ls[6]].split(":")), ls[6]
                    )
                )
    # annovar annotate
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -geneanno "
        f"-dbtype refGene -hgvs -buildver hg38 --outfile {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_geneanno.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype gnomad_exome --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_gnomad.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype gnomad_exome --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_gnomad.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype clinvar --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_clinvar.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype cosmicTier12 --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_cosmic.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype intervar --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_intervar.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"{dic['annovar']['path']}/annotate_variation.pl -filter "
        f"-dbtype oncokb --thread {dic['samtools']['threads']} "
        f"-buildver hg38 -out {outd}/Mid-files/09 "
        f"{outd}/Mid-files/08_avinput_left.xls {dic['annovar']['humandb']} "
        f"2> {outd}/logs/09_oncokb.stderr"
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"python3 {dic['scripts']}/merge_annotations.py "
        f"-oncokb_tsv {dic['docs']}/cancerGeneList.tsv -hotspots  "
        f"{dic['docs']}/hotspots_twoColumns.bed -outdir {outd}/Mid-files "
        f"-fq-prefix {sample} "
    )
    # print(cmd)
    os.system(cmd)
    os.system(
        f"rm -f {outd}/Mid-files/06_trimPrimer_* {outd}/Mid-files/07_variant_* "
        f"{outd}/Mid-files/*_filtered {outd}/Mid-files/09.log "
    )
    # filter variant
    if dic["sample-type"] == "Tissue":
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"python3 {dic['scripts']}/variant_QC_filter.py "
            f"-snvIndels {outd}/{sample}_SNV_Indels_total.xls -depth 300 "
            f"-alt 5 -frequency 1 -out {outd}/{sample}_SNV_Indels_filter.xls "
        )
        # print(cmd)
        os.system(cmd)
    else:
        cmd = (
            f"singularity exec -B {mount_paras} {dic['image']} "
            f"python3 {dic['scripts']}/variant_QC_filter.py "
            f"-snvIndels {outd}/{sample}_SNV_Indels_total.xls -depth 500 "
            f"-alt 5 -frequency 0.5 -out {outd}/{sample}_SNV_Indels_filter.xls "
        )
        # print(cmd)
        os.system(cmd)
    filter_df = pd.read_csv(
        f"{outd}/{sample}_SNV_Indels_filter.xls", sep="\t", index_col=None
    )
    filter_df = filter_df.sort_values(by="Frequency(%)", ascending=False)
    filter_df.to_csv(f"{outd}/{sample}_SNV_Indels_filter.xls", sep="\t", index=None)


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
    return methy_dic, amplicon_dic


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
        (
            cpg
            / (
                cpg
                + v["CT-CT"]["CpG_unmethylated"]
                + v["GA-GA"]["CpG_unmethylated"]
                + v["CT-GA"]["CpG_unmethylated"]
                + v["GA-CT"]["CpG_unmethylated"]
            )
            * 100
        )
        if (
            cpg
            + v["CT-CT"]["CpG_unmethylated"]
            + v["GA-GA"]["CpG_unmethylated"]
            + v["CT-GA"]["CpG_unmethylated"]
            + v["GA-CT"]["CpG_unmethylated"]
        )
        != 0
        else 0
    )
    methylated["nonCpG_methylated%"] = (
        (
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
        )
        if (
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
        != 0
        else 0
    )
    methylated["methylated_all%"] = (
        (
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
        )
        if (
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
        != 0
        else 0
    )
    # CT GA strand
    cpg_ct = v["CT-CT"]["CpG_methylated"] + v["GA-CT"]["CpG_methylated"]
    cpg_ga = v["GA-GA"]["CpG_methylated"] + v["CT-GA"]["CpG_methylated"]
    chg_ct = v["CT-CT"]["CHG_methylated"] + v["GA-CT"]["CHG_methylated"]
    chg_ga = v["GA-GA"]["CHG_methylated"] + v["CT-GA"]["CHG_methylated"]
    chh_ct = v["CT-CT"]["CHH_methylated"] + v["GA-CT"]["CHH_methylated"]
    chh_ga = v["GA-GA"]["CHH_methylated"] + v["CT-GA"]["CHH_methylated"]
    methylated["CpG_methylated_CT%"] = (
        (
            cpg_ct
            / (cpg_ct + v["CT-CT"]["CpG_unmethylated"] + v["GA-CT"]["CpG_unmethylated"])
            * 100
        )
        if (cpg_ct + v["CT-CT"]["CpG_unmethylated"] + v["GA-CT"]["CpG_unmethylated"])
        != 0
        else 0
    )
    methylated["CpG_methylated_GA%"] = (
        (
            cpg_ga
            / (cpg_ga + v["CT-GA"]["CpG_unmethylated"] + v["GA-GA"]["CpG_unmethylated"])
            * 100
        )
        if (cpg_ga + v["CT-GA"]["CpG_unmethylated"] + v["GA-GA"]["CpG_unmethylated"])
        != 0
        else 0
    )
    methylated["nonCpG_methylated_CT%"] = (
        (
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
        )
        if (
            chg_ct
            + chh_ct
            + v["GA-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHG_unmethylated"]
            + v["CT-CT"]["CHH_unmethylated"]
            + v["GA-CT"]["CHH_unmethylated"]
        )
        != 0
        else 0
    )
    methylated["nonCpG_methylated_GA%"] = (
        (
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
        )
        if (
            chg_ga
            + chh_ga
            + v["GA-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHG_unmethylated"]
            + v["CT-GA"]["CHH_unmethylated"]
            + v["GA-GA"]["CHH_unmethylated"]
        )
        != 0
        else 0
    )
    methylated["Methylation_CT%"] = (
        (
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
        )
        if (
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
        != 0
        else 0
    )
    methylated["Methylation_GA%"] = (
        (
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
        )
        if (
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
        != 0
        else 0
    )
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


def summary_methylate(outd, dic, sample, mount_paras):
    # fq aspect
    stat_dic = {}
    try:
        qc_df = pd.read_json(f"{outd}/Mid-files/01_CleanFq.json")
    except:
        print("Wrong: Can`t find Sample-Dir/Mid-files/01_CleanFq.json. ")
        exit()
    raw_reads_num = int(qc_df["read1_before_filtering"]["total_reads"])
    clean_reads_num = int(qc_df["read1_after_filtering"]["total_reads"])
    raw_G = round(raw_reads_num * 300 / 1000000000, 2)
    if raw_G > 3:
        print(
            f"Warning: this sample size is too big ({raw_G} G), "
            f"may take a long time to calculate the depth. Please be patient and wait"
        )
    short_read = int(qc_df["filtering_result"]["too_short_reads"])
    total_read = int(qc_df["summary"]["before_filtering"]["total_reads"])
    short_rate = round(float(short_read / total_read) * 100, 2)
    clean_rate = round(clean_reads_num / raw_reads_num * 100, 2)
    q30 = qc_df["summary"]["before_filtering"]["q30_rate"]
    # bwa aspect
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools flagstat -@ 2 {outd}/Mid-files/03_bismark_sort.bam >{outd}/Mid-files/zz.stat.xls1"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools flagstat -@ 2 {outd}/Mid-files/05_dedup_sort.bam >{outd}/Mid-files/zz.stat.xls2"
    )
    std = subprocess.run(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"""awk '$0~"properly paired"{{print $1}}' {outd}/Mid-files/zz.stat.xls1""",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    mapping_reads_num = int(std.stdout.strip()) / 2
    map_rate = round(mapping_reads_num / clean_reads_num * 100, 2)
    std = subprocess.run(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"""awk '$0~"properly paired"{{print $1}}' {outd}/Mid-files/zz.stat.xls2""",
        shell=True,
        capture_output=True,
        encoding="utf-8",
    )
    unique_read = int(int(std.stdout.strip()) / 2)
    unique_percend = round(unique_read / mapping_reads_num * 100, 2)
    # primer dedup
    bed_dic, methy_dic, amplicon_dic = {}, {}, {}
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
    with open(dic["bed"], "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            bed_dic[ke] = {}
            bed_dic[ke]["chr"] = ls[0]
            bed_dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            bed_dic[ke]["primer"] = ls[3]
            bed_dic[ke]["length"] = int(ls[4])
            bed_dic[ke]["strand"] = ls[5]
            bed_dic[ke]["orient"] = "CT" if ls[5] == "+" else "GA"
            bed_dic[ke]["depth"] = 0
            amplicon_dic[ls[3]] = {
                "cpg_unmeth": 0,
                "cpg_meth": 0,
                "noncpg_unmeth": 0,
                "noncpg_meth": 0,
            }
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -@ 2 -b -f 64 -F 2048 "
        f"{outd}/Mid-files/05_dedup_sort.bam >{outd}/Mid-files/zz.stat.bam"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/zz.stat.bam"
    )
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/Mid-files/zz.stat.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in bed_dic.keys():
                    if (
                        bed_dic[k]["chr"] == chr_name
                        and bed_dic[k]["strand"] == read_strand
                        and bed_dic[k]["orient"] == read.get_tag("XG")
                        and np.abs(position - bed_dic[k]["pos"]) <= 2
                        and match_len > bed_dic[k]["length"]
                    ):
                        primername = bed_dic[k]["primer"]
                        bed_dic[k]["depth"] += 1
                        fit_read += 1
                        convert = read.get_tag("XR") + "-" + read.get_tag("XG")
                        methy_dic, amplicon_dic = GetMethylationDict(
                            methy_dic,
                            convert,
                            read.get_tag("XM"),
                            primername,
                            amplicon_dic,
                        )
    uni_on_target_read = round(fit_read / unique_read * 100, 2)
    depth_dic, fit = {}, 0
    for k in bed_dic.keys():
        if bed_dic[k]["primer"] not in depth_dic.keys():
            depth_dic[bed_dic[k]["primer"]] = {}
            depth_dic[bed_dic[k]["primer"]]["Unique depth(X)"] = bed_dic[k]["depth"]
        else:
            depth_dic[bed_dic[k]["primer"]]["Unique depth(X)"] += bed_dic[k]["depth"]
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
    depth_5dedup_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    # sort3 bam
    bed_dic = {}
    with open(dic["bed"], "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            ke = f"{ls[0]}:{ls[1]}:{ls[2]}"
            bed_dic[ke] = {}
            bed_dic[ke]["chr"] = ls[0]
            bed_dic[ke]["pos"] = int(ls[1]) if ls[5] == "+" else int(ls[2])
            bed_dic[ke]["primer"] = ls[3]
            bed_dic[ke]["length"] = int(ls[4])
            bed_dic[ke]["strand"] = ls[5]
            bed_dic[ke]["orient"] = "CT" if ls[5] == "+" else "GA"
            bed_dic[ke]["depth"] = 0
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -@ 2 -b -f 64 -F 2048 "
        f"{outd}/Mid-files/03_bismark_sort.bam >{outd}/Mid-files/zz.stat.bam"
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/zz.stat.bam"
    )
    fit_read = 0
    with pysam.AlignmentFile(f"{outd}/Mid-files/zz.stat.bam", "rb") as f:
        for read in f.fetch():
            cigartuples = read.cigartuples
            read_strand = "-" if read.is_reverse else "+"
            match_len = count_match_length(read_strand, cigartuples)  # method
            if match_len > 0:
                chr_name = read.reference_name
                chr_sta = int(read.reference_start)
                position = chr_sta + match_len if read.is_reverse else chr_sta
                for k in bed_dic.keys():
                    if (
                        bed_dic[k]["chr"] == chr_name
                        and bed_dic[k]["strand"] == read_strand
                        and bed_dic[k]["orient"] == read.get_tag("XG")
                        and np.abs(position - bed_dic[k]["pos"]) <= 2
                        and match_len > bed_dic[k]["length"]
                    ):
                        bed_dic[k]["depth"] += 1
                        fit_read += 1
    on_target_rate = round(fit_read / mapping_reads_num * 100, 2)
    depth_dic, fit = {}, 0
    for k in bed_dic.keys():
        if bed_dic[k]["primer"] not in depth_dic.keys():
            depth_dic[bed_dic[k]["primer"]] = {}
            depth_dic[bed_dic[k]["primer"]]["Depth(X)"] = bed_dic[k]["depth"]
        else:
            depth_dic[bed_dic[k]["primer"]]["Depth(X)"] += bed_dic[k]["depth"]
    mean_depth = int(fit_read / len(depth_dic))
    for k in depth_dic.keys():
        if depth_dic[k]["Depth(X)"] >= mean_depth * 0.2:
            fit += 1
    uniformity = round(fit / len(depth_dic) * 100, 2)
    depth_3sort_df = (
        pd.DataFrame(depth_dic).T.reset_index().rename(columns={"index": "Primer"})
    )
    depth_df = pd.merge(depth_3sort_df, depth_5dedup_df, on="Primer", how="outer")
    depth_df["Unique depth(X)"] = depth_df["Unique depth(X)"].astype(int)
    depth_df.loc[
        depth_df["Depth(X)"] == 0,
        ["CpG converted%", "NonCpG converted%", "Conversion all%"],
    ] = "/"
    depth_df.to_csv(f"{outd}/{sample}_depth.xls", sep="\t", index=None)

    # final aspect
    stat_dic["Sample"] = sample
    stat_dic["Raw(G)"] = raw_G
    stat_dic["Raw read(PE)"] = raw_reads_num
    stat_dic["<50 bp read%"] = short_rate
    stat_dic["Clean read%"] = clean_rate
    stat_dic["Q30 rate%"] = round(float(q30) * 100, 2)
    stat_dic["Mapping rate%"] = map_rate
    stat_dic["Unique read(PE)"] = unique_read
    stat_dic["Unique percent%"] = unique_percend
    stat_dic["Unique on target rate%"] = uni_on_target_read
    stat_dic["On target rate%"] = on_target_rate
    stat_dic["Unique mean depth(X)"] = uniquemeandepth
    stat_dic["Mean depth(X)"] = mean_depth
    stat_dic["Unique uniformity%"] = UniqueUniformity
    stat_dic["Uniformity%"] = uniformity
    stat_dic = GetMethylationRate(stat_dic, methy_dic)  # method
    stat_df = pd.DataFrame([stat_dic])
    stat_df.to_csv(f"{outd}/{sample}_basic_summary.xls", sep="\t", index=None)


def methylate_1to5(outd, dic, sample, mount_paras):
    # fastp
    fq1 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R1_suffix"]}'
    fq2 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R2_suffix"]}'
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"fastp -i {fq1} -I {fq2} -o {outd}/Mid-files/01_CleanFq_1.fq.gz "
        f"-O {outd}/Mid-files/01_CleanFq_2.fq.gz "
        f"--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        f"--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--correction --fix_mgi_id --cut_tail -l 30 -w {dic['fastp']['threads']} "
        f"--reads_to_process {dic['fastp']['read_to_process']} "
        f"-j {outd}/Mid-files/01_CleanFq.json --overlap_len_require 10 "
        f"-h {outd}/Mid-files/01_CleanFq.html 2> {outd}/logs/01_fastp.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # add umi
    counta = 0
    with gzip.open(
        f"{outd}/Mid-files/02_UmiFq_1.fq.gz", "wt", compresslevel=4
    ) as fq1_output:
        with gzip.open(
            f"{outd}/Mid-files/02_UmiFq_2.fq.gz", "wt", compresslevel=4
        ) as fq2_output:
            with gzip.open(f"{outd}/Mid-files/01_CleanFq_1.fq.gz", "rt") as handle1:
                with gzip.open(f"{outd}/Mid-files/01_CleanFq_2.fq.gz", "rt") as handle2:
                    for tu in zip(handle1, handle2):
                        counta += 1
                        if counta % 4 == 1:
                            title1 = tu[0].split()[0]
                            title2 = tu[1].split()[0]
                        if counta % 4 == 2:
                            tmp_umi = ":" + tu[1][:20]
                            fq1_output.write(title1 + tmp_umi + "\n")
                            fq1_output.write(tu[0])
                            fq2_output.write(title2 + tmp_umi + "\n")
                            fq2_output.write(tu[1])
                        if counta % 4 == 3:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
                        if counta % 4 == 0:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
    # bismark
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bismark -1 {outd}/Mid-files/02_UmiFq_1.fq.gz -2 {outd}/Mid-files/02_UmiFq_2.fq.gz "
        f"-B {outd}/Mid-files/03_bismark --genome {dic['bismark']['genome']} "
        f"--temp_dir {outd}/Mid-files --non_directional --local -p {dic['bismark']['threads']} "
        f"1> {outd}/logs/03_bismark.stdout 2> {outd}/logs/03_bismark.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # deduplicate_bismark
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"deduplicate_bismark {outd}/Mid-files/03_bismark_pe.bam "
        f"-o 05 --output_dir {outd}/Mid-files --barcode  "
        f"2> {outd}/logs/05_deduplicate.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # samtools sort
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -bh -f 3 -q 10 -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/03_bismark_pe.bam "
        f"-o {outd}/Mid-files/03_bismark_pe2.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools sort -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/03_bismark_pe2.bam -o {outd}/Mid-files/03_bismark_sort.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/03_bismark_sort.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools view -bh -f 3 -q 10 -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/05.deduplicated.bam "
        f"-o {outd}/Mid-files/05.deduplicated2.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools sort -@ {dic['samtools']['threads']} "
        f"{outd}/Mid-files/05.deduplicated2.bam -o {outd}/Mid-files/05_dedup_sort.bam "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index {outd}/Mid-files/05_dedup_sort.bam "
    )
    # print(cmd)
    os.system(cmd)
    summary_methylate(outd, dic, sample, mount_paras)  # method
    os.system(f"rm -f {outd}/Mid-files/03_bismark_pe* {outd}/Mid-files/05.dedup* ")


def fusion_read_summary(outd, dic, sample, mount_paras):
    # depth count
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bedtools intersect -a {outd}/Mid-files/03_sort_filter.bam -b {dic['bed']} "
        f"-F 0.9 -bed -wa -wb -s > {outd}/Mid-files/zz.3.bam.intersect "
    )
    os.system(
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"bedtools intersect -a {outd}/Mid-files/07_dedup_filter.bam -b {dic['bed']} "
        f"-F 0.9 -bed -wa -wb -s > {outd}/Mid-files/zz.7.bam.intersect "
    )
    # dispose bed file
    primer_info = {}
    with open(dic["bed"], "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            primer_info[ls[3]] = ls
    # primer depth
    read_primer = {}
    with open(dic["bed"], "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            read_primer[ls[3]] = {}
            read_primer[ls[3]]["Depth"] = 0
            read_primer[ls[3]]["Unique depth"] = 0
    with open(f"{outd}/Mid-files/zz.3.bam.intersect", "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            mapped_start = int(ls[1])
            mapped_end = int(ls[2])
            primer_start = int(ls[13])
            primer_end = int(ls[14])
            primer_name = ls[15]
            mapped_strand = ls[5]
            read_pair = int(ls[3].split("/")[1])
            if mapped_strand == "+":
                if read_pair == 1:
                    if abs(mapped_start - primer_start) <= 2:
                        read_primer[primer_name]["Depth"] += 1
                else:
                    if abs(mapped_end - primer_end) <= 2:
                        read_primer[primer_name]["Depth"] += 1
            else:
                if read_pair == 1:
                    if abs(mapped_end - primer_end) <= 2:
                        read_primer[primer_name]["Depth"] += 1
                else:
                    if abs(mapped_start - primer_start) <= 2:
                        read_primer[primer_name]["Depth"] += 1
    with open(f"{outd}/Mid-files/zz.7.bam.intersect", "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            mapped_start = int(ls[1])
            mapped_end = int(ls[2])
            primer_start = int(ls[13])
            primer_end = int(ls[14])
            primer_name = ls[15]
            mapped_strand = ls[5]
            if mapped_strand == "+":
                if abs(mapped_start - primer_start) <= 2:
                    read_primer[primer_name]["Unique depth"] += 1
            else:
                if abs(mapped_end - primer_end) <= 2:
                    read_primer[primer_name]["Unique depth"] += 1
    depth_df = (
        pd.DataFrame(read_primer).T.reset_index().rename(columns={"index": "Primer"})
    )
    depth_df.to_csv(f"{outd}/{sample}_depth.xls", sep="\t", index=None)
    # fq aspect
    stat_dic = {}
    try:
        qc_df = pd.read_json(f"{outd}/Mid-files/01_CleanFq.json")
    except:
        print("Wrong: Can`t find Sample-Dir/Mid-files/01_CleanFq.json ")
        exit()
    raw_reads_num = int(qc_df["read1_before_filtering"]["total_reads"])
    clean_reads_num = int(qc_df["read1_after_filtering"]["total_reads"])
    raw_G = round(raw_reads_num * 300 / 1000000000, 2)
    short_read = int(qc_df["filtering_result"]["too_short_reads"])
    total_read = int(qc_df["summary"]["before_filtering"]["total_reads"])
    short_rate = round(float(short_read / total_read) * 100, 2)
    clean_rate = round(clean_reads_num / raw_reads_num * 100, 2)
    q30 = qc_df["summary"]["before_filtering"]["q30_rate"]
    # alignment aspect
    with open(f"{outd}/Mid-files/03_star.Log.final.out", "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            if "Number of input reads" in line:
                total_reads = int(ls[1])
            elif "Uniquely mapped reads number" in line:
                unique_mapped_reads = int(ls[1])
            elif "Number of reads mapped to multiple loci" in line:
                multi_loci_reads = int(ls[1])
            elif "Number of reads mapped to too many loci" in line:
                many_loci_reads = int(ls[1])
            elif "Number of chimeric reads" in line:
                chimeric_reads = int(ls[1])
    map_read = unique_mapped_reads + multi_loci_reads + many_loci_reads + chimeric_reads
    map_rate = round(map_read / clean_reads_num * 100, 2)
    with open(f"{outd}/Mid-files/07_star.Log.final.out", "r") as handle:
        for line in handle:
            ls = line.strip().split("\t")
            if "Number of input reads" in line:
                unique_total_reads = int(ls[1])
            elif "Uniquely mapped reads number" in line:
                unique_mapped_reads = int(ls[1])
            elif "Number of reads mapped to multiple loci" in line:
                multi_loci_reads = int(ls[1])
            elif "Number of reads mapped to too many loci" in line:
                many_loci_reads = int(ls[1])
            elif "Number of chimeric reads" in line:
                chimeric_reads = int(ls[1])
    uni_map_read = (
        unique_mapped_reads + multi_loci_reads + many_loci_reads + chimeric_reads
    )
    uni_percent = round(unique_total_reads / total_reads * 100, 2)
    # stat
    control = ["CHMP2A", "GPI", "RAB7A", "VCP"]
    uni_mean_depth = int(depth_df["Unique depth"].sum() / len(depth_df))
    mean_depth = int(depth_df["Depth"].sum() / len(depth_df))
    uni_control_mean_depth = int(
        depth_df[depth_df["Primer"].str.split("_").str[0].isin(control)][
            "Unique depth"
        ].sum()
        / len(
            depth_df[depth_df["Primer"].str.split("_").str[0].isin(control)][
                "Unique depth"
            ]
        )
    )
    uni_on_target_rate = round(depth_df["Unique depth"].sum() / uni_map_read * 100, 2)
    on_target_rate = round(depth_df["Depth"].sum() / map_read * 100, 2)
    uni_uniformity = round(
        len(depth_df[depth_df["Unique depth"] > uni_mean_depth * 0.2])
        / len(depth_df)
        * 100,
        2,
    )
    uniformity = round(
        len(depth_df[depth_df["Depth"] > mean_depth * 0.2]) / len(depth_df) * 100, 2
    )
    # final aspect
    stat_dic["Sample"] = sample
    stat_dic["Raw(G)"] = raw_G
    stat_dic["Raw read(PE)"] = raw_reads_num
    stat_dic["<50 bp read%"] = short_rate
    stat_dic["Clean read%"] = clean_rate
    stat_dic["Q30 rate%"] = round(float(q30) * 100, 2)
    stat_dic["Mapping rate%"] = map_rate
    stat_dic["Unique read(PE)"] = uni_map_read
    stat_dic["Unique percent%"] = uni_percent
    stat_dic["Unique on target rate%"] = uni_on_target_rate
    stat_dic["On target rate%"] = on_target_rate
    stat_dic["Unique mean depth(X)"] = uni_mean_depth
    stat_dic["Mean depth(X)"] = mean_depth
    stat_dic["Unique uniformity%"] = uni_uniformity
    stat_dic["Uniformity%"] = uniformity
    stat_dic["Control mean depth"] = uni_control_mean_depth
    stat_df = pd.DataFrame([stat_dic])
    stat_df.to_csv(f"{outd}/{sample}_basic_summary.xls", sep="\t", index=None)


def fusion_1to5(outd, dic, sample, mount_paras):
    # fastp
    fq1 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R1_suffix"]}'
    fq2 = f'{dic["samples"][sample]["directory"]}/{sample}{dic["samples"][sample]["R2_suffix"]}'
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"fastp -i {fq1} -I {fq2} -o {outd}/Mid-files/01_CleanFq_1.fq.gz "
        f"-O {outd}/Mid-files/01_CleanFq_2.fq.gz "
        f"--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        f"--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--correction --fix_mgi_id --cut_tail -l 50 -w {dic['fastp']['threads']} "
        f"--umi --umi_loc read2 --umi_len 9 --umi_prefix UMI "
        f"--reads_to_process {dic['fastp']['read_to_process']} "
        f"-j {outd}/Mid-files/01_CleanFq.json --overlap_len_require 10 "
        f"-h {outd}/Mid-files/01_CleanFq.html 2> {outd}/logs/01_fastp.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # add umi
    counta = 0
    with gzip.open(
        f"{outd}/Mid-files/02_UmiFq_1.fq.gz", "wt", compresslevel=4
    ) as fq1_output:
        with gzip.open(
            f"{outd}/Mid-files/02_UmiFq_2.fq.gz", "wt", compresslevel=4
        ) as fq2_output:
            with gzip.open(f"{outd}/Mid-files/01_CleanFq_1.fq.gz", "rt") as handle1:
                with gzip.open(f"{outd}/Mid-files/01_CleanFq_2.fq.gz", "rt") as handle2:
                    for tu in zip(handle1, handle2):
                        counta += 1
                        if counta % 4 == 1:
                            title1 = tu[0].split()[0]
                            title2 = tu[1].split()[0]
                        if counta % 4 == 2:
                            tmp_umi = tu[1][:11]
                            fq1_output.write(title1 + tmp_umi + "\n")
                            fq1_output.write(tu[0])
                            fq2_output.write(title2 + tmp_umi + "\n")
                            fq2_output.write(tu[1])
                        if counta % 4 == 3:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
                        if counta % 4 == 0:
                            fq1_output.write(tu[0])
                            fq2_output.write(tu[1])
    # load fusion genome
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR --genomeLoad LoadAndExit --genomeDir {dic['STAR']['genomeDir']} "
        f"--outFileNamePrefix {outd}/.genomeLoad/ 1> /dev/null "
    )
    # print(cmd)
    os.system(cmd)
    # STAR software alignment
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR --readFilesIn {outd}/Mid-files/02_UmiFq_1.fq.gz {outd}/Mid-files/02_UmiFq_2.fq.gz "
        f"--runThreadN {dic['STAR']['threads']} --genomeDir {dic['STAR']['genomeDir']} "
        f"--readFilesCommand zcat --outSAMstrandField intronMotif --outSAMunmapped Within "
        f"--alignInsertionFlush Right --outSAMtype BAM SortedByCoordinate --outSAMattrRGline  "
        f' "ID:flowcell.lane"  --genomeLoad LoadAndKeep --chimOutType Junctions WithinBAM SoftClip '
        f"--quantMode GeneCounts --chimSegmentMin 12 --chimJunctionOverhangMin 8 "
        f"--chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 "
        f"--alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 "
        f"--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --chimScoreJunctionNonGTAG -4 "
        f"--chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 "
        f"--limitBAMsortRAM 10000000000 --outFileNamePrefix {outd}/Mid-files/03_star. "
        f"1> {outd}/logs/03_star.stdout 2> {outd}/logs/03_star.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # samtools index
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index -@ {dic['samtools']['threads']} {outd}/Mid-files/03_star.Aligned.sortedByCoord.out.bam "
        f"1> {outd}/logs/03_samtools_index.stdout"
    )
    # print(cmd)
    os.system(cmd)
    # gencore
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"gencore -r {dic['bwa']['fasta']} -u UMI -s 1 -i {outd}/Mid-files/03_star.Aligned.sortedByCoord.out.bam "
        f"-o {outd}/Mid-files/05_dedup_not_sort.bam -h {outd}/Mid-files/05_gencore.html "
        f"-j {outd}/Mid-files/05_gencore.json 2> {outd}/logs/05_gencore.stderr"
    )
    # print(cmd)
    os.system(cmd)
    # samtools sort
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools sort -@ {dic['samtools']['threads']} {outd}/Mid-files/05_dedup_not_sort.bam "
        f"> {outd}/Mid-files/05_dedup_sort.bam "
    )
    # print(cmd)
    os.system(cmd)
    # samtools index
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index -@ {dic['samtools']['threads']} {outd}/Mid-files/05_dedup_sort.bam "
        f"1> {outd}/logs/05_samtools_index.stdout"
    )
    # print(cmd)
    os.system(cmd)
    # bamtofastq
    read1_sequence, read1_quality, read2_sequence, read2_quality, check_pair = (
        {},
        {},
        {},
        {},
        {},
    )
    with pysam.AlignmentFile(f"{outd}/Mid-files/05_dedup_sort.bam", "rb") as handle:
        for read in handle:
            quality = "".join(map(lambda x: chr(x + 33), read.get_forward_qualities()))
            if read.qname not in check_pair:
                check_pair[read.qname] = {}
            if read.is_read1:
                check_pair[read.qname]["read1"] = ""
                read1_sequence[read.qname] = read.get_forward_sequence()
                read1_quality[read.qname] = quality
            elif read.is_read2:
                check_pair[read.qname]["read2"] = ""
                read2_sequence[read.qname] = read.get_forward_sequence()
                read2_quality[read.qname] = quality
    with gzip.open(f"{outd}/Mid-files/06_pairend_1.fq.gz", "wt") as read1:
        with gzip.open(f"{outd}/Mid-files/06_pairend_2.fq.gz", "wt") as read2:
            for read in check_pair:
                if len(check_pair[read]) == 2:
                    read1.write(
                        "@{0}\n{1}\n+\n{2}\n".format(
                            read, read1_sequence[read], read1_quality[read]
                        )
                    )
                    read2.write(
                        "@{0}\n{1}\n+\n{2}\n".format(
                            read, read2_sequence[read], read2_quality[read]
                        )
                    )
    # single end fastp
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"fastp -i {outd}/Mid-files/06_pairend_1.fq.gz -I {outd}/Mid-files/06_pairend_2.fq.gz "
        f"--merge --merged_out {outd}/Mid-files/06_merge_se.fq.gz "
        f"--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        f"--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT "
        f"--correction -l 50 -w {dic['fastp']['threads']} "
        f"-j {outd}/Mid-files/06_merge.json --overlap_len_require 10 "
        f"-h {outd}/Mid-files/06_merge.html 2> {outd}/logs/06_fastp.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # STAR software alignment twice
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR --readFilesIn {outd}/Mid-files/06_merge_se.fq.gz "
        f"--runThreadN {dic['STAR']['threads']} --genomeDir {dic['STAR']['genomeDir']} "
        f"--readFilesCommand zcat --outSAMstrandField intronMotif --outSAMunmapped Within "
        f"--alignInsertionFlush Right --outSAMtype BAM SortedByCoordinate --outSAMattrRGline "
        f' "ID:flowcell.lane"  --genomeLoad LoadAndKeep --chimOutType Junctions WithinBAM SoftClip '
        f"--quantMode GeneCounts --chimSegmentMin 12 --chimJunctionOverhangMin 8 "
        f"--chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 "
        f"--alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 "
        f"--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --chimScoreJunctionNonGTAG -4 "
        f"--chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 "
        f"--limitBAMsortRAM 10000000000 --outFileNamePrefix {outd}/Mid-files/07_star. --outSAMattributes All "
        f"1> {outd}/logs/07_star.stdout 2> {outd}/logs/07_star.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # samtools index
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"samtools index -@ {dic['samtools']['threads']} {outd}/Mid-files/07_star.Aligned.sortedByCoord.out.bam "
        f"2> {outd}/logs/07_samtools_index.stderr"
    )
    # print(cmd)
    os.system(cmd)
    # remove fusion genome
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR --genomeLoad Remove --genomeDir {dic['STAR']['genomeDir']} "
        f"--outFileNamePrefix {outd}/.genomeLoad/ 1> /dev/null "
    )
    # print(cmd)
    os.system(cmd)
    # filter bam
    with pysam.AlignmentFile(
        f"{outd}/Mid-files/03_star.Aligned.sortedByCoord.out.bam", "rb"
    ) as handle:
        with pysam.AlignmentFile(
            f"{outd}/Mid-files/03_sort_filter.bam", "wb", template=handle
        ) as output:
            for read in handle:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.is_read2 and not read.is_supplementary:
                    continue
                if read.is_read1:
                    if read.is_reverse:
                        if read.cigartuples[-1][0] == 4:
                            if read.cigartuples[-1][1] > 2:
                                continue
                        output.write(read)
                    else:
                        if read.cigartuples[0][0] == 4:
                            if read.cigartuples[0][1] > 2:
                                continue
                        output.write(read)
                else:
                    if read.is_reverse:
                        if read.cigartuples[0][0] == 4:
                            if read.cigartuples[0][1] > 2:
                                continue
                        output.write(read)
                    else:
                        if read.cigartuples[-1][0] == 4:
                            if read.cigartuples[-1][1] > 2:
                                continue
                        output.write(read)
    with pysam.AlignmentFile(
        f"{outd}/Mid-files/07_star.Aligned.sortedByCoord.out.bam", "rb"
    ) as handle:
        with pysam.AlignmentFile(
            f"{outd}/Mid-files/07_dedup_filter.bam", "wb", template=handle
        ) as output:
            for read in handle:
                if read.is_unmapped or read.is_secondary:
                    continue
                if read.is_reverse:
                    if read.cigartuples[-1][0] == 4:
                        if read.cigartuples[-1][1] > 2:
                            continue
                    output.write(read)
                else:
                    if read.cigartuples[0][0] == 4:
                        if read.cigartuples[0][1] > 2:
                            continue
                    output.write(read)
    # read summary
    fusion_read_summary(outd, dic, sample, mount_paras)  # method
    # STAR-Fusion
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR-Fusion --genome_lib_dir {dic['star_fusion']['star_fusion_lib']} "
        f"--chimeric_junction {outd}/Mid-files/07_star.Chimeric.out.junction --output_dir {outd}/Mid-files "
        f"--examine_coding_effect --min_sum_frags 1 --CPU {dic['star_fusion']['threads']} "
        f"2> {outd}/logs/08_star_fusion.stderr"
    )
    # print(cmd)
    os.system(cmd)
    # STAR_to_cancer_introns
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"STAR_to_cancer_introns.py --ctat_genome_lib {dic['star_fusion']['star_fusion_lib']} "
        f"--chimJ_file {outd}/Mid-files/07_star.Chimeric.out.junction "
        f"--bam_file {outd}/Mid-files/07_star.Aligned.sortedByCoord.out.bam "
        f"--SJ_tab_file {outd}/Mid-files/07_star.SJ.out.tab --output_prefix {outd}/Mid-files/08_ctat "
        f"2> {outd}/logs/08_star_intron.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # filter threshold
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"python3 {dic['scripts']}/fusion_QC_filter.py "
        f"-fusion_db {dic['star_fusion']['star_fusion_lib']} -predict_tsv "
        f"{outd}/Mid-files/star-fusion.fusion_predictions.tsv -oncokb_tsv "
        f"{dic['docs']}/cancerGeneList.tsv -output_bedpe {outd}/Mid-files/09_on_target.bedpe "
        f"-outdir {outd}/Mid-files -prefix_name {sample} "
        f"2> {outd}/logs/09_qc_filter.stderr "
    )
    # print(cmd)
    os.system(cmd)
    # igv
    cmd = (
        f"singularity exec -B {mount_paras} {dic['image']} "
        f"create_report {outd}/Mid-files/09_on_target.bedpe "
        f"{dic['star_fusion']['star_fusion_lib']}/ref_genome.fa "
        f"--track-config {outd}/Mid-files/09_trackConfigs.json "
        f"--output {outd}/{sample}_fusion_igv.html "
        f"2> {outd}/logs/10_igv.stderr "
    )
    # print(cmd)
    os.system(cmd)


def make_cnv_baseline(outd, bam_str, dic, mount_paras, image):
    # bait
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"cnvkit.py target {outd}/zz.region.bed --annotate {dic['docs']}/refFlat.txt "
        f"-o {outd}/zz.tmp.bait 2> {outd}/zz.log"
    )
    # print(cmd)
    os.system(cmd)
    # bait bed
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"python3 {dic['scripts']}/mark_bed_cnvfusion.py --input {outd}/zz.tmp.bait --cnv-gene "
        f"{dic['docs']}/general_CNVgene.txt --output {outd}/zz.tmp.bait.bed"
    )
    # print(cmd)
    os.system(cmd)
    # cnv baseline making
    baselinename = "CNV" + dic["baseline-name"]
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"cnvkit.py batch -n {bam_str} -m amplicon -f {dic['bwa']['fasta']} "
        f"-t {outd}/zz.tmp.bait.bed --target-avg-size 50 "
        f"-d {outd}/BaselineCreator-mid-files/CNVcalling --output-reference {outd}/{baselinename} "
        f"2> {outd}/zz.log"
    )
    # print(cmd)
    os.system(cmd)


def make_msi_baseline(outd, bam_str, dic, mount_paras, image, sample_dic):
    with open(f"{outd}/zz.configure.xls", "w") as fi:
        for fq in sample_dic.keys():
            fi.write(
                "{}\t{}\n".format(
                    fq,
                    f"{outd}/BaselineCreator-mid-files/{fq}/Mid-files/{fq}_dedup_sort.bam",
                )
            )
    # get bed microsatellites.tsv
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"""awk -F "\\t" -v OFS="\\t" 'NR>1{{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' """
        f"{dic['docs']}/GRCh38.microsatellites.tsv >{outd}/zz.msi.tsv.nohead1"
    )
    # print(cmd)
    os.system(cmd)
    # intersect
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"bedtools intersect -a {outd}/zz.msi.tsv.nohead1 -b "
        f"{outd}/zz.region.bed -wa |cut -f 1,3- >{outd}/zz.msi.tsv.nohead"
    )
    # print(cmd)
    os.system(cmd)
    #
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"head -n 1 {dic['docs']}/GRCh38.microsatellites.tsv|cat - {outd}/zz.msi.tsv.nohead >{outd}/zz.msi.tsv"
    )
    # print(cmd)
    os.system(cmd)
    # msi baseline making
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"msisensor-pro baseline -d {outd}/zz.msi.tsv -i {outd}/zz.configure.xls "
        f"-o {outd}/BaselineCreator-mid-files/MSIcalling -c 15 1> {outd}/zz.log"
    )
    # print(cmd)
    os.system(cmd)
    # designate name
    baselinename = "MSI" + dic["baseline-name"]
    cmd = (
        f"singularity exec -B {mount_paras} {image} "
        f"cp {outd}/BaselineCreator-mid-files/MSIcalling/zz.msi.tsv_baseline {outd}/{baselinename}"
    )
    # print(cmd)
    os.system(cmd)


def ParameterCombine(
    sample_dic,
    outd,
    variantCall,
    cnvCall,
    msiCall,
    methCall,
    fusionCall,
    mount_paras,
    dic,
    image,
    autoRemove,
    baselineCall,
    baselineType,
):
    primer_convert(outd, dic, mount_paras)  # method
    if variantCall or cnvCall or msiCall:
        print("Run mutation analysis ... ...")
        for sample in sample_dic.keys():
            sampledir = Path(f"{outd}/{sample}/logs")
            mid_dir = Path(f"{outd}/{sample}/Mid-files")
            if not sampledir.exists():
                sampledir.mkdir(parents=True)
                mid_dir.mkdir(parents=True)
            variant_1to5(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            # read summary
            summary_variant(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            if cnvCall:
                print("CNV calling ... ...")
                variant_cnv(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            if msiCall:
                print("MSI calling ... ...")
                variant_msi(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            if variantCall:
                print("Mutation Annotating  ... ...")
                variant_6to11(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            os.system(f"rm -f {outd}/{sample}/zz* {outd}/{sample}/Mid-files/zz*")

            with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                yaml.dump(
                    dic,
                    fi,
                    default_flow_style=False,
                    allow_unicode=True,
                    sort_keys=False,
                )
            if autoRemove:
                midfile = Path(f"{outd}/{sample}/Mid-files")
                if midfile.exists():
                    os.system(f"rm -rf {outd}/{sample}/Mid-files")
        os.system(f"rm -f {outd}/AtomSeq_config.yaml {outd}/zz* ")
    if methCall:
        print("Run Atomseq target methylation ... ...")
        for sample in sample_dic.keys():
            sampledir = Path(f"{outd}/{sample}/logs")
            mid_dir = Path(f"{outd}/{sample}/Mid-files")
            if not sampledir.exists():
                sampledir.mkdir(parents=True)
                mid_dir.mkdir(parents=True)
            methylate_1to5(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            os.system(f"rm -f {outd}/{sample}/zz* {outd}/{sample}/Mid-files/zz*")

            with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                yaml.dump(
                    dic,
                    fi,
                    default_flow_style=False,
                    allow_unicode=True,
                    sort_keys=False,
                )
            if autoRemove:
                midfile = Path(f"{outd}/{sample}/Mid-files")
                if midfile.exists():
                    os.system(f"rm -rf {outd}/{sample}/Mid-files")
        os.system(f"rm -f {outd}/AtomSeq_config.yaml {outd}/zz* ")
    if fusionCall:
        print("Run fusion RNA data ... ...")
        for sample in sample_dic.keys():
            sampledir = Path(f"{outd}/{sample}/logs")
            mid_dir = Path(f"{outd}/{sample}/Mid-files")
            if not sampledir.exists():
                sampledir.mkdir(parents=True)
                mid_dir.mkdir(parents=True)
            fusion_1to5(f"{outd}/{sample}", dic, sample, mount_paras)  # method
            os.system(
                f"rm -rf {outd}/{sample}/.genomeLoad/ {outd}/{sample}/zz* {outd}/{sample}/Mid-files/zz*"
            )

            with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                yaml.dump(
                    dic,
                    fi,
                    default_flow_style=False,
                    allow_unicode=True,
                    sort_keys=False,
                )
            if autoRemove:
                midfile = Path(f"{outd}/{sample}/Mid-files")
                if midfile.exists():
                    os.system(f"rm -rf {outd}/{sample}/Mid-files")
        os.system(f"rm -f {outd}/AtomSeq_config.yaml {outd}/zz* ")
    if baselineCall:
        print("Creating new baseline file ... ...")
        bam_str = ""
        for sample in sample_dic.keys():
            sampledir = Path(f"{outd}/BaselineCreator-mid-files/{sample}/logs")
            mid_dir = Path(f"{outd}/BaselineCreator-mid-files/{sample}/Mid-files")
            if not sampledir.exists():
                sampledir.mkdir(parents=True)
                mid_dir.mkdir(parents=True)
            variant_1to5(
                f"{outd}/BaselineCreator-mid-files/{sample}", dic, sample, mount_paras
            )  # method
            bam_str += f"{outd}/BaselineCreator-mid-files/{sample}/Mid-files/{sample}_dedup_sort.bam "
        if baselineType == "CNV":
            make_cnv_baseline(outd, bam_str, dic, mount_paras, image)  # method
        else:
            make_msi_baseline(
                outd, bam_str, dic, mount_paras, image, sample_dic
            )  # method

        with open(f"{outd}/Sample_config.yaml", "w") as fi:
            yaml.dump(
                dic,
                fi,
                default_flow_style=False,
                allow_unicode=True,
                sort_keys=False,
            )
        if autoRemove:
            midfile = Path(f"{outd}/BaselineCreator-mid-files")
            if midfile.exists():
                os.system(f"rm -rf {outd}/BaselineCreator-mid-files")
        os.system(f"rm -f {outd}/AtomSeq_config.yaml {outd}/zz* ")


def defaultYaml(
    sample_dic,
    config_dic,
    variantCall,
    cnvCall,
    msiCall,
    methCall,
    fusionCall,
    sampleType,
    cnvBaseline,
    msiBaseline,
    primerbed,
    supportRead,
    outd,
    image,
    autoRemove,
    dataSize,
    threads,
    baselineCall,
    baselineType,
    baselineName,
):
    if dataSize:
        read_num = int(dataSize * 1000000000 / (2 * 150) + 1)
    else:
        read_num = int(dataSize)
    varopt = True if variantCall else False
    cnvopt = True if cnvCall else False
    msiopt = True if msiCall else False
    methopt = True if methCall else False
    fusopt = True if fusionCall else False
    baselineCreatoropt = True if baselineCall else False
    sampleCharater = sampleType if sampleType else False
    cnvBline = cnvBaseline if cnvBaseline else False
    msiBline = msiBaseline if msiBaseline else False
    baselineType = baselineType if baselineType else False
    fasta_dir = os.path.dirname(f"{config_dic['Hg38_Fasta_Path']}")
    default = {
        "samples": sample_dic,
        "python3": "python3",
        "scripts": "/snakemake/bin/dev",
        "docs": "/snakemake/resources",
        # "scripts": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v2.8.1/bin/dev",
        # "docs": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v2.8.1/resources",
        "bed": primerbed,
        "image": image,
        "variant-call": varopt,
        "cnv-call": cnvopt,
        "msi-call": msiopt,
        "cnv-baseline": cnvBline,
        "msi-baseline": msiBline,
        "methylation-call": methopt,
        "fusion-call": fusopt,
        "sample-type": sampleCharater,
        "baseline-creator": baselineCreatoropt,
        "baseline-output-type": baselineType,
        "baseline-name": baselineName,
        "auto-remove": autoRemove,
        "bedtools": "bedtools",
        "fastp": {
            "path": "fastp",
            "read_to_process": read_num,
            "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "threads": threads,
        },
        "bwa": {
            "path": "bwa",
            "fasta": config_dic["Hg38_Fasta_Path"],
            "bwa_index": f"{fasta_dir}/Bwa_index/hg38",
            # "bwa_index": config_dic["Bwa_Index_Path"],
            "read_group": r' "@RG\tID:flowcell1.lane1\tLB:library1\tSM:sample\tPL:ILLUMINA" ',
            "threads": threads,
        },
        "bismark": {
            "path": "bismark",
            "genome": fasta_dir,
            # "genome": config_dic["Bismark_Index_Path"],
            "threads": threads,
        },
        "STAR": {
            "path": "STAR",
            "genomeDir": f"{config_dic['Fusion_data_library']}/ref_genome.fa.star.idx",
            "threads": threads,
        },
        "star_fusion": {
            "path": "STAR-Fusion",
            "star_fusion_lib": config_dic["Fusion_data_library"],
            "threads": threads,
        },
        "deduplicate_bismark": {"path": "deduplicate_bismark"},
        "samtools": {"path": "samtools", "threads": threads},
        "gatk": {
            "path": "gatk",
            "known_snps": f"{config_dic['Variant_library']}/GATK_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "known_indels": f"{config_dic['Variant_library']}/GATK_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        },
        "gencore": {"path": "gencore", "support_reads": supportRead},
        "cnvkit": {
            "path": "cnvkit.py",
            "Rscript_path": "/snakemake/bin/dev/Rscript",
            # "Rscript_path": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/bin/dev/Rscript",
        },
        "msisensor_pro": {"path": "msisensor-pro"},
        "varscan2": {"path": "varscan", "min_coverage": 40, "min_ave_quality": 11},
        "annovar": {
            "path": f"{config_dic['Variant_library']}/annovar",
            "humandb": f"{config_dic['Variant_library']}/annovar_library",
            "threads": threads,
        },
        "vep": {
            "path": "vep",
            "vep_dir": "/opt/conda/envs/vep/share/ensembl-vep-108.2-0",
            # "vep_dir": "/biotron/Application/Software/install/mambaforge/envs/vep108/share/ensembl-vep-108.2-0",
            "vep_library": f"{config_dic['Variant_library']}/VEP_library",
        },
        "snpEff": {
            "path": f"{config_dic['Variant_library']}/snpEff/snpEff.jar",
        },
    }
    with open(f"{outd}/AtomSeq_config.yaml", "w") as fi:
        yaml.dump(
            default, fi, default_flow_style=False, allow_unicode=True, sort_keys=False
        )
    return default


def func(
    image,
    fqdir,
    fq_lis,
    outd,
    configFile,
    primerbed,
    variantCall,
    sampleType,
    cnvCall,
    cnvBaseline,
    msiCall,
    msiBaseline,
    methCall,
    fusionCall,
    baselineCall,
    baselineType,
    baselineName,
    threads,
    dataSize,
    supportRead,
    autoRemove,
):
    pydir = os.path.dirname(os.path.realpath(__file__))
    if not image:
        image = glob.glob(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "*.sif")
        )[0]
    else:
        image = os.path.abspath(image)

    # check abspath
    softlink, hardlink = [], []
    for fq in glob.glob(f"{fqdir}/*_*1.f*q.gz"):
        sample = fq.split("/")[-1].split("_")[0]
        if os.path.islink(fq):
            softlink_path = os.path.dirname(fq)
            original_path = os.path.dirname(os.readlink(fq))
            if sample in fq_lis:
                if softlink_path not in softlink:
                    softlink.append(softlink_path)
                if original_path not in hardlink:
                    hardlink.append(original_path)
        else:
            original_path = os.path.dirname(fq)
            if sample in fq_lis:
                if original_path not in hardlink:
                    hardlink.append(original_path)
    mount_path = hardlink.copy()
    outd = os.path.abspath(outd)
    configFile = os.path.abspath(configFile)
    primerbed = os.path.abspath(primerbed)

    outputdir = Path(outd)
    if not outputdir.exists():
        outputdir.mkdir(parents=True)

    mount_path += [outd, configFile, primerbed]

    if cnvBaseline:
        cnvBaseline = os.path.abspath(cnvBaseline)
        mount_path.append(cnvBaseline)
    if msiBaseline:
        msiBaseline = os.path.abspath(msiBaseline)
        mount_path.append(msiBaseline)
    if softlink:
        mount_path += softlink
    mount_paras = " -B ".join(mount_path)

    config_dic, mount_paras = GetConfigFile(
        configFile, mount_paras, image, pydir
    )  # method
    sample_dic = FastqDispose(fqdir, fq_lis)  # method
    CheckModeConflict(
        variantCall,
        cnvCall,
        msiCall,
        methCall,
        fusionCall,
        sampleType,
        msiBaseline,
        baselineCall,
        baselineType,
    )  # method
    defaultDic = defaultYaml(
        sample_dic,
        config_dic,
        variantCall,
        cnvCall,
        msiCall,
        methCall,
        fusionCall,
        sampleType,
        cnvBaseline,
        msiBaseline,
        primerbed,
        supportRead,
        outd,
        image,
        autoRemove,
        dataSize,
        threads,
        baselineCall,
        baselineType,
        baselineName,
    )  # method
    ParameterCombine(
        sample_dic,
        outd,
        variantCall,
        cnvCall,
        msiCall,
        methCall,
        fusionCall,
        mount_paras,
        defaultDic,
        image,
        autoRemove,
        baselineCall,
        baselineType,
    )  # method

    time_end = datetime.datetime.now()
    print(
        f"Total time spent: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    )


def main():
    argv = argparse_line()
    func(
        argv["image"],
        argv["fq_dir"],
        argv["fq_prefix"],
        argv["outdir"],
        argv["config"],
        argv["bed"],
        argv["variant"],
        argv["sample_type"],
        argv["cnv"],
        argv["cnv_baseline"],
        argv["msi"],
        argv["msi_baseline"],
        argv["meth"],
        argv["fusion"],
        argv["baselineCreator"],
        argv["baseline_type"],
        argv["baseline_name"],
        argv["threads"],
        argv["data_to_analyse"],
        argv["supporting_reads"],
        argv["auto_remove"],
    )


if __name__ == "__main__":
    main()
