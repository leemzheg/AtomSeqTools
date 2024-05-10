#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/03/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os
import argparse, textwrap
import datetime
from pathlib import Path

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="Make hg38 index for variant and methylation",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-image", metavar="STR", help="Singularity image [*.sif]", required=True
    )
    required.add_argument(
        "-fasta", metavar="STR", help="Input hg38.fasta", required=True
    )
    argv = vars(parser.parse_args())
    return argv


def func(image, fasta):
    # pydir = os.path.dirname(os.path.realpath(__file__))
    outd = os.path.dirname(fasta)
    cmd = (
        f"singularity exec -B {outd} {image} "
        f"samtools dict {fasta} >{outd}/hg38.dict"
    )
    # print(cmd)
    os.system(cmd)
    bwa_dir = Path(f"{outd}/Bwa_index")
    if not bwa_dir.exists():
        bwa_dir.mkdir(parents=True)
    cmd = (
        f"singularity exec -B {outd} {image} "
        f"bwa index {fasta} -p {outd}/Bwa_index/hg38 "
    )
    # print(cmd)
    os.system(cmd)
    cmd = (
        f"singularity exec -B {outd} {image} "
        f"bismark_genome_preparation --path_to_aligner /opt/conda/bin --verbose {outd}"
    )
    # print(cmd)
    os.system(cmd)

    time_end = datetime.datetime.now()
    print(
        f"Make hg38 index used: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    )


def main():
    argv = argparse_line()
    func(
        argv["image"],
        argv["fasta"],
    )


if __name__ == "__main__":
    main()
