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
    required.add_argument("-oncokb_tsv", metavar="i", help="i", required=True)
    required.add_argument("-hotspots", metavar="i", help="i", required=True)
    required.add_argument("-outdir", metavar="o", help="o", required=True)
    required.add_argument("-fq-prefix", metavar="i", help="i", required=True)
    argv = vars(parser.parse_args())
    return argv


def read_oncokb(tsv):
    oncokb_dic = {}
    with open(tsv) as f:
        for line in f:
            ls = line.strip().split("\t")
            oncokb_dic[ls[4]] = ""
            key = ls[5].split(".")[0]
            oncokb_dic[key] = ""
    return oncokb_dic


def left_annotate(hotspot, outdir, oncokb_dic):
    right_dic = {}
    hotspot_dic = {}
    varient_lis = [
        "exonic",
        "exonic;splicing",
        "ncRNA_exonic",
        "ncRNA_exonic;splicing",
        "splicing",
    ]
    with open(hotspot) as f:
        for line in f:
            ls = line.strip().split("\t")
            key = f"{ls[0]}:{ls[1]}"
            hotspot_dic[key] = ""
    # exonic variants
    with open(f"{outdir}/09.exonic_variant_function") as f:
        for line in f:
            ls = line.strip().split("\t")
            right_dic[ls[-2]] = {}
            col3_lis = ls[2].split(",")
            if f"{ls[3]}:{ls[4]}" in hotspot_dic:
                right_dic[ls[-2]]["hotspot"] = "HOT"
            else:
                right_dic[ls[-2]]["hotspot"] = "Non"
            for item in col3_lis[:-1]:
                NM = item.split(":")
                if len(NM) < 5:
                    continue
                if NM[1] in oncokb_dic.keys():
                    right_dic[ls[-2]][
                        "geneanno"
                    ] = f"{NM[0]}:{ls[1]}:{NM[2]}:{NM[1]}:{NM[3]}:{NM[4]}"
            if "geneanno" not in right_dic[ls[-2]].keys():
                right_dic[ls[-2]]["geneanno"] = "Non:Non:Non:Non:Non:Non"
    # ncRNA and splicing variants
    with open(f"{outdir}/09.variant_function") as f:
        for line in f:
            ls = line.strip().split("\t")
            if ls[-2] not in right_dic.keys() and ls[0] in varient_lis:
                right_dic[ls[-2]] = {}
                right_dic[ls[-2]]["geneanno"] = "Non:Non:Non:Non:Non:Non"
                if f"{ls[2]}:{ls[3]}" in hotspot_dic.keys():
                    right_dic[ls[-2]]["hotspot"] = "HOT"
                else:
                    right_dic[ls[-2]]["hotspot"] = "Non"
    # clinvar cosmic gnomAD intervar anotate
    database_lis = ["clinvar", "cosmicTier12", "gnomad_exome", "intervar", "oncokb"]
    for db in database_lis:
        tmp_dic = {}
        with open(f"{outdir}/09.hg38_{db}_dropped") as f:
            for line in f:
                ls = line.strip().split("\t")
                tmp_dic[ls[-2]] = ls[1]
        for k in right_dic.keys():
            if k in tmp_dic.keys():
                right_dic[k][db] = tmp_dic[k]
            else:
                right_dic[k][db] = "Non"

    return right_dic


def right_annotate(outdir, oncokb_dic, right_dic):
    # oncokb mafannotator dictory
    # with open(f"{outdir}/8_target_oncokb.maf") as f:
    #     for line in f:
    #         ls = line.strip().split("\t")
    #         if re.match(r"Hugo_", line):
    #             continue
    #         if ls[102] in right_dic.keys():
    #             right_dic[ls[102]]["oncokb"] = ls[117] + ":" + ls[119]
    # snpEff dictory
    with open(f"{outdir}/07_variant.snpeff.vcf") as f:
        for line in f:
            ls = line.strip().split("\t")
            if re.match(r"#", line) or ls[2] not in right_dic.keys():
                continue
            col8_lis = ls[7].split(";")[5].split(",")
            for item in col8_lis:
                NM = item.split("|")
                if (
                    NM[6].split(".")[0] in oncokb_dic.keys()
                    and "hgvs" not in right_dic[ls[2]]
                ):
                    if len(NM[10].split("p.")[-1]) == 0:
                        suffix = "?"
                    else:
                        suffix = NM[10].split("p.")[-1]
                    right_dic[ls[2]][
                        "hgvs"
                    ] = f'{NM[3]}:{NM[1]}:exon{NM[8].split("/")[0]}:{NM[6].split(".")[0]}:{NM[9]}:p.{suffix}'
            if "hgvs" not in right_dic[ls[2]]:
                NM = col8_lis[0].split("|")
                if len(NM[10].split("p.")[-1]) == 0:
                    suffix = "?"
                else:
                    suffix = NM[10].split("p.")[-1]
                right_dic[ls[2]][
                    "hgvs"
                ] = f'{NM[3]}:{NM[1]}:exon{NM[8].split("/")[0]}:{NM[6].split(".")[0]}:{NM[9]}:p.{suffix}'

    return right_dic


def func(oncokb, hotspot, outdir, samplename):
    # oncokb representative and annotate
    oncokb_dic = read_oncokb(oncokb)
    right_dic = left_annotate(hotspot, outdir, oncokb_dic)
    final_dic = right_annotate(outdir, oncokb_dic, right_dic)
    title = (
        "Chrom\tStart\tEnd\tRef\tAlt\tTotal_count\tAlt_count"
        "\tFrequency(%)\tMutation_type\tAnnotation\tHGVS.p"
        "\tONCOGENIC:MUTATION_EFFECT\tGnomAD_exome_EAS\tCLNALLELEID:CLNSIG"
        "\tInterVar\tCosmicTier12\tHotspots\n"
    )
    # write into file
    with open(f"{outdir}/08_avinput_fix.xls") as f, open(
        f"{outdir}/../{samplename}_SNV_Indels_total.xls", "w"
    ) as fi:
        fi.write(f"{title}")
        for line in f:
            ls = line.strip().split("\t")
            if ls[-1] in final_dic.keys():
                if final_dic[ls[-1]]["geneanno"] == "Non:Non:Non:Non:Non:Non":
                    anno = final_dic[ls[-1]]["hgvs"].split(":")
                    if re.search(r"splic", anno[1]):
                        anno[1] = "splicing mutation"
                    annotation = anno[:1] + anno[2:]
                    fi.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            "\t".join(ls[:-1]),
                            anno[1],
                            ":".join(annotation),
                            anno[-1],
                            final_dic[ls[-1]]["oncokb"],
                            final_dic[ls[-1]]["gnomad_exome"],
                            final_dic[ls[-1]]["clinvar"],
                            final_dic[ls[-1]]["intervar"],
                            final_dic[ls[-1]]["cosmicTier12"],
                            final_dic[ls[-1]]["hotspot"],
                        )
                    )
                else:
                    anno = final_dic[ls[-1]]["geneanno"].split(":")
                    hgvs_c = final_dic[ls[-1]]["hgvs"].split(":")[-2]
                    hgvs_p = final_dic[ls[-1]]["hgvs"].split(":")[-1]
                    if re.search(r"splic", anno[1]):
                        anno[1] = "splicing mutation"
                    annotation = f"{anno[0]}:{anno[2]}:{anno[3]}:{hgvs_c}:{anno[-1]}"
                    fi.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            "\t".join(ls[:-1]),
                            anno[1],
                            annotation,
                            hgvs_p,
                            final_dic[ls[-1]]["oncokb"],
                            final_dic[ls[-1]]["gnomad_exome"],
                            final_dic[ls[-1]]["clinvar"],
                            final_dic[ls[-1]]["intervar"],
                            final_dic[ls[-1]]["cosmicTier12"],
                            final_dic[ls[-1]]["hotspot"],
                        )
                    )


def main():
    argv = argparse_line()
    func(
        argv["oncokb_tsv"],
        argv["hotspots"],
        argv["outdir"],
        argv["fq_prefix"],
    )


if __name__ == "__main__":
    main()
