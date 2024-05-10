#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/03/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os, yaml
import argparse, textwrap
import datetime, glob
from pathlib import Path

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
        help="    Specify fastq prefix to be analyzed(can be many, delimited by space)\n"
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
        help="    If call variant, it`s necessary to select sample type. Possible values: {cfDNA, Tissue}",
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
    meth = parser.add_argument_group("Fusion options")
    meth.add_argument(
        "-fusion",
        action="store_true",
        help="    Enable pipeline to analyse RNA fusion data",
    )
    optional = parser.add_argument_group("Other options")
    optional.add_argument(
        "-threads",
        metavar="INT",
        help="    Number of threads to use [default:12]",
        type=int,
        default=12,
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


def GetConfigFile(config, mount_paras, image):
    config_dic = {}
    pydir = os.path.dirname(os.path.realpath(__file__))
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
    variantCall, cnvCall, msiCall, methCall, fusionCall, sampleType, msiBaseline
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
                "Methylation mode at the same time(Only one mode can be activated), please check"
            )
            exit()
        if fusionCall:
            print(
                "Traceback: You have already opened Variant/cnv/msi mode, you can`t open "
                "Fusion mode at the same time(Only one mode can be activated), please check"
            )
            exit()
    if methCall:
        if fusionCall:
            print(
                "Traceback: You have already opened Methylation mode, you can`t open "
                "Fusion mode at the same time(Only one mode can be activated), please check"
            )
            exit()
        if variantCall or cnvCall or msiCall:
            print(
                "Traceback: You have already opened Methylation mode, you can`t open "
                "Variant/cnv/msi mode at the same time(Only one mode can be activated), please check"
            )
            exit()
    if fusionCall:
        if methCall:
            print(
                "Traceback: You have already opened Fusion mode, you can`t open "
                "Methylation mode at the same time(Only one mode can be activated), please check"
            )
            exit()
        if variantCall or cnvCall or msiCall:
            print(
                "Traceback: You have already opened Methylation mode, you can`t open "
                "Variant/cnv/msi mode at the same time(Only one mode can be activated), please check"
            )
            exit()


def ParameterCombine(
    sample_dic,
    config_dic,
    dataSize,
    outd,
    primerbed,
    variantCall,
    cnvCall,
    msiCall,
    methCall,
    fusionCall,
    threads,
    supportRead,
    mount_paras,
    fqdir,
    fq_lis,
    defaultDic,
    image,
    autoRemove,
):
    if variantCall or cnvCall or msiCall:
        print("Run variant preprocessing ... ...")
        if variantCall:
            cmd = (
                f"singularity exec -B {mount_paras} {image} "
                f"snakemake -s /snakemake/bin/snakefile_variant "
                f"--configfile {outd}/AtomSeq_config.yaml --cores {threads} --rerun-incomplete "
                f"-d {outd} -q progress --use-conda --printshellcmds"
            )
            # print(f"   {cmd}")
            os.system(cmd)
            for sample in sample_dic.keys():
                with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                    yaml.dump(
                        defaultDic,
                        fi,
                        default_flow_style=False,
                        allow_unicode=True,
                        sort_keys=False,
                    )
                    if autoRemove:
                        midfile = Path(f"{outd}/{sample}/Mid-files")
                        if midfile.exists():
                            os.system(f"rm -rf {outd}/{sample}/Mid-files")
            os.system(f"rm -rf {outd}/.snakemake " f"{outd}/AtomSeq_config.yaml ")
        elif cnvCall:
            cmd = (
                f"singularity exec -B {mount_paras} {image} "
                f"snakemake -s /snakemake/bin/snakefile_cnv "
                f"--configfile {outd}/AtomSeq_config.yaml --cores {threads} --rerun-incomplete "
                f"-d {outd} -q progress --use-conda --printshellcmds"
            )
            # print(f"   {cmd}")
            os.system(cmd)
            for sample in sample_dic.keys():
                with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                    yaml.dump(
                        defaultDic,
                        fi,
                        default_flow_style=False,
                        allow_unicode=True,
                        sort_keys=False,
                    )
                    if autoRemove:
                        midfile = Path(f"{outd}/{sample}/Mid-files")
                        if midfile.exists():
                            os.system(f"rm -rf {outd}/{sample}/Mid-files")
            os.system(f"rm -rf {outd}/.snakemake " f"{outd}/AtomSeq_config.yaml ")
        elif msiCall:
            cmd = (
                f"singularity exec -B {mount_paras} {image} "
                f"snakemake -s /snakemake/bin/snakefile_msi "
                f"--configfile {outd}/AtomSeq_config.yaml --cores {threads} --rerun-incomplete "
                f"-d {outd} -q progress --use-conda --printshellcmds"
            )
            # print(f"   {cmd}")
            os.system(cmd)
            for sample in sample_dic.keys():
                with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                    yaml.dump(
                        defaultDic,
                        fi,
                        default_flow_style=False,
                        allow_unicode=True,
                        sort_keys=False,
                    )
                    if autoRemove:
                        midfile = Path(f"{outd}/{sample}/Mid-files")
                        if midfile.exists():
                            os.system(f"rm -rf {outd}/{sample}/Mid-files")
            os.system(f"rm -rf {outd}/.snakemake " f"{outd}/AtomSeq_config.yaml ")
    if methCall:
        print("Run atomseq target methylation ... ...")
        cmd = (
            f"singularity exec -B {mount_paras} {image} "
            f"snakemake -s /snakemake/bin/snakefile_methylation "
            f"--configfile {outd}/AtomSeq_config.yaml --cores {threads} --rerun-incomplete "
            f"-d {outd} -q progress --use-conda --printshellcmds"
        )
        os.system(cmd)
        for sample in sample_dic.keys():
            with open(f"{outd}/{sample}/Sample_config.yaml", "w") as fi:
                yaml.dump(
                    defaultDic,
                    fi,
                    default_flow_style=False,
                    allow_unicode=True,
                    sort_keys=False,
                )
                if autoRemove:
                    midfile = Path(f"{outd}/{sample}/Mid-files")
                    if midfile.exists():
                        os.system(f"rm -rf {outd}/{sample}/Mid-files")
        os.system(f"rm -rf {outd}/.snakemake " f"{outd}/AtomSeq_config.yaml ")
    if fusionCall:
        print("Run RNA fusion data ... ...")
        if autoRemove:
            cmd = (
                f"singularity exec -B {mount_paras} {image} "
                f"/snakemake/bin/prereq_atomseq-rna-fusion.py --image {image} "
                f"--fastq-dir {fqdir} --call-variants-lib-dir {config_dic['Fusion_data_library']} "
                f"--primer-bed {primerbed} --output-dir {outd} --min-junction-reads 5 "
                f"--min-fusion-frequency 0.01 --min-primer-depth 50 --supporting-reads {supportRead} "
                f"--fastq-prefix {' '.join(fq_lis)} --data-size {dataSize} --auto-remove"
            )
        else:
            cmd = (
                f"singularity exec -B {mount_paras} {image} "
                f"/snakemake/bin/prereq_atomseq-rna-fusion.py --image {image} "
                f"--fastq-dir {fqdir} --call-variants-lib-dir {config_dic['Fusion_data_library']} "
                f"--primer-bed {primerbed} --output-dir {outd} --min-junction-reads 5 "
                f"--min-fusion-frequency 0.01 --min-primer-depth 50 --supporting-reads {supportRead} "
                f"--fastq-prefix {' '.join(fq_lis)} --data-size {dataSize}"
            )
        # print(f"   {cmd}")
        os.system(cmd)


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
    sampleCharater = sampleType if sampleType else False
    cnvBline = cnvBaseline if cnvBaseline else False
    msiBline = msiBaseline if msiBaseline else False
    fasta_dir = os.path.dirname(f"{config_dic['Hg38_Fasta_Path']}")
    default = {
        "samples": sample_dic,
        "python3": "python3",
        "scripts": "/snakemake/bin/dev",
        "docs": "/snakemake/resources",
        # "scripts": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/bin/dev",
        # "docs": "/biotron/Application/WorkStation/limengzheng/atomSeqTools/AtomSeqTools_v1.0.8/resources",
        "bed": primerbed,
        "image": image,
        "variant-call": varopt,
        "cnv-call": cnvopt,
        "msi-call": msiopt,
        "methylation-call": methopt,
        "fusion-call": fusopt,
        "sample-type": sampleCharater,
        "cnv-baseline": cnvBline,
        "msi-baseline": msiBline,
        "auto-remove": autoRemove,
        "bedtools": "bedtools",
        "fastp": {
            "path": "fastp",
            "read_to_process": read_num,
            "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "threads": 8,
        },
        "bwa": {
            "path": "bwa",
            "fasta": config_dic["Hg38_Fasta_Path"],
            "bwa_index": f"{fasta_dir}/Bwa_index/hg38",
            # "bwa_index": config_dic["Bwa_Index_Path"],
            "read_group": r' "@RG\tID:flowcell1.lane1\tLB:library1\tSM:sample\tPL:ILLUMINA" ',
            "threads": 8,
        },
        "bismark": {
            "path": "bismark",
            "genome": fasta_dir,
            # "genome": config_dic["Bismark_Index_Path"],
            "threads": 8,
        },
        "deduplicate_bismark": {"path": "deduplicate_bismark"},
        "samtools": {"path": "samtools", "threads": 8},
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
            "threads": 8,
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
    threads,
    dataSize,
    supportRead,
    autoRemove,
):
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

    config_dic, mount_paras = GetConfigFile(configFile, mount_paras, image)  # method
    sample_dic = FastqDispose(fqdir, fq_lis)  # method
    CheckModeConflict(
        variantCall, cnvCall, msiCall, methCall, fusionCall, sampleType, msiBaseline
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
    )  # method
    ParameterCombine(
        sample_dic,
        config_dic,
        dataSize,
        outd,
        primerbed,
        variantCall,
        cnvCall,
        msiCall,
        methCall,
        fusionCall,
        threads,
        supportRead,
        mount_paras,
        fqdir,
        fq_lis,
        defaultDic,
        image,
        autoRemove,
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
        argv["threads"],
        argv["data_to_analyse"],
        argv["supporting_reads"],
        argv["auto_remove"],
    )


if __name__ == "__main__":
    main()
