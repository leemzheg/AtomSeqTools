#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# @Time : 2024/03/01 9:00
# @Author : Li Mengzheng
# @E-mail : mengzheng-li@ebiotron.com

import re, os, json
import argparse, textwrap
import pandas as pd
import datetime
import gffutils, pysam

time_start = datetime.datetime.now()


def argparse_line():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawTextHelpFormatter
    )
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "-predict_tsv",
        metavar="i",
        help="input star fusion prediction tsv file",
        required=True,
    )
    required.add_argument(
        "-fusion_db", metavar="i", help="input fusion database", required=True
    )
    required.add_argument(
        "-oncokb_tsv", metavar="i", help="input oncokb tsv", required=True
    )
    required.add_argument(
        "-outdir", metavar="o", help="output directory", required=True
    )
    required.add_argument(
        "-output_bedpe", metavar="o", help="output bam file", required=True
    )
    required.add_argument(
        "-prefix_name", metavar="str", help="sample name", required=True
    )
    # optional = parser.add_argument_group("Optional arguments")
    # optional.add_argument(
    #     "-threads",
    #     metavar="INT",
    #     help="Number of threads to use [default:12]",
    #     type=int,
    #     default=12,
    # )
    argv = vars(parser.parse_args())
    return argv


def read_ctat_splicing(ctat_cancer_introns, ctat_introns):
    target = {"METx14del": "MET_Exon14_Skipping", "EGFRvIII": "EGFR_Exon2-7_Skipping"}
    exon_skipping = {}
    exon_range = {}
    with open(ctat_cancer_introns, "r") as handle:
        for line in handle:
            newline = line.rstrip().split("\t")
            if re.match("intron", line):
                continue
            variant_name = newline[7]
            read_num = int(newline[3])
            loc = newline[0]
            chromosome = loc.split(":")[0]
            start = int(loc.split(":")[-1].split("-")[0])
            end = int(loc.split(":")[-1].split("-")[-1])
            if variant_name in target:
                exon_skipping[target[variant_name]] = {}
                exon_skipping[target[variant_name]]["readnum"] = read_num
                exon_skipping[target[variant_name]]["loc"] = loc
                exon_skipping[target[variant_name]]["totalread"] = 0
                exon_range[target[variant_name]] = {}
                exon_range[target[variant_name]][chromosome] = [start, end]
    with open(ctat_introns, "r") as handle:
        for line in handle:
            newline = line.rstrip().split("\t")
            if re.match("intron", line):
                continue
            chromosome = newline[0].split(":")[0]
            start = int(newline[0].split(":")[-1].split("-")[0])
            end = int(newline[0].split(":")[-1].split("-")[-1])
            uniq_mapped = int(newline[3])
            for skipping in exon_range:
                if chromosome in exon_range[skipping]:
                    if (
                        exon_range[skipping][chromosome][0] <= start
                        and exon_range[skipping][chromosome][1] >= end
                    ):
                        exon_skipping[skipping]["totalread"] += uniq_mapped
    return exon_skipping


def match_readid_primer(outd):
    read_primer, readid_primer_dic = {}, {}
    with open(f"{outd}/zz.7.bam.intersect", "r") as handle:
        for line in handle:
            newline = line.strip().split("\t")
            mapped_start = int(newline[1])
            mapped_end = int(newline[2])
            read_id = newline[3]
            primer_start = int(newline[13])
            primer_end = int(newline[14])
            primer_name = newline[15]
            mapped_strand = newline[5]
            if mapped_strand == "+":
                if abs(mapped_start - primer_start) <= 2:
                    if read_id not in read_primer:
                        read_primer[read_id] = {}
                        read_primer[read_id][primer_name] = ""
                    else:
                        read_primer[read_id][primer_name] = ""
            else:
                if abs(mapped_end - primer_end) <= 2:
                    if read_id not in read_primer:
                        read_primer[read_id] = {}
                        read_primer[read_id][primer_name] = ""
                    else:
                        read_primer[read_id][primer_name] = ""
    for read in read_primer:
        for primer in read_primer[read]:
            readid_primer_dic[read] = primer

    return readid_primer_dic


def calculate_primer_frequency(junction_reads, read_primer):
    fusion_primer_depth = {}
    fusion_primer_read = {}
    fusion_intarget_read = []
    fusion_outoftarget_read = []
    for read in junction_reads.split(","):
        read_name = read.split("@")[-1]
        if read_name in read_primer:
            fusion_intarget_read.append(read_name)
            primer = read_primer[read_name]
            if primer not in fusion_primer_depth:
                fusion_primer_depth[primer] = 1
                fusion_primer_read[primer] = [read_name]
            else:
                fusion_primer_depth[primer] += 1
                fusion_primer_read[primer].append(read_name)
        else:
            fusion_outoftarget_read.append(read_name)

    return (
        fusion_primer_depth,
        fusion_primer_read,
        fusion_intarget_read,
        fusion_outoftarget_read,
    )


def find_exon_number(break_point, sqldb, oncokb_accession):
    kit_accession = []
    db = gffutils.FeatureDB(sqldb, keep_order=True)

    break_point = break_point.split(":")
    region = "{0}:{1}-{1}".format(break_point[0], break_point[1])
    strand = break_point[2]

    candidate_transcript = {
        "Kit_Select": "",
        "OncoKB_Select": "",
        "MANE_Select": "",
        "RefSeq_Select": "",
    }
    for i in db.region(region, strand=strand, featuretype=["exon", "intron"]):
        transcript = "".join(i.attributes["transcript_id"]).split(".")[0]
        if i.featuretype == "exon":
            number = "".join(i.attributes["exon_number"])
        else:
            number = i.attributes["exon_number"][0]
        #####################################################################
        if transcript in kit_accession:
            candidate_transcript["Kit_Select"] = (
                transcript + ":" + i.featuretype + number
            )
        if transcript in oncokb_accession:
            candidate_transcript["OncoKB_Select"] = (
                transcript + ":" + i.featuretype + number
            )
        if "tag" in i.attributes:
            if "MANE Select" in i.attributes["tag"]:
                candidate_transcript["MANE_Select"] = (
                    transcript + ":" + i.featuretype + number
                )
            if "RefSeq Select" in i.attributes["tag"]:
                candidate_transcript["RefSeq_Select"] = (
                    transcript + ":" + i.featuretype + number
                )

    exon_number = []
    for select in candidate_transcript:
        if candidate_transcript[select]:
            exon_number.append(select + ":" + candidate_transcript[select])

    if not exon_number:
        primary_exon_num = "NA"
        exon_number = "NA"
    else:
        primary_exon_num = exon_number[0].split(":")[-1]
        exon_number = ",".join(exon_number)

    return primary_exon_num, exon_number


def extract_fusion_intarget_read_align(
    input_bam_file,
    output_bam_file,
    all_fusion_intarget_read,
    ncbi_txt,
    igv_json_file,
    skipping,
):

    all_read_id = []
    with pysam.AlignmentFile(input_bam_file, "rb") as handle:
        with pysam.AlignmentFile(output_bam_file, "wb", template=handle) as output:
            for gene in skipping:
                skipping_reads = handle.fetch(
                    skipping[gene]["chr"],
                    skipping[gene]["start"],
                    skipping[gene]["end"],
                )
                for read in skipping_reads:
                    if read not in all_read_id:
                        output.write(read)
                        all_read_id.append(read)
            for read in handle:
                if read.query_name in all_fusion_intarget_read:
                    if read not in all_read_id:
                        output.write(read)
                        all_read_id.append(read)

    sorted_bam = "{0}.sorted.bam".format(output_bam_file.split(".bam")[0])
    pysam.sort("-o", sorted_bam, output_bam_file)
    pysam.index(sorted_bam)

    track_config = [
        {
            "name": "Refseq Genes",
            "format": "refgene",
            "url": ncbi_txt,
        },
        {
            "name": "Alignments",
            "url": "{0}".format(sorted_bam),
            "indexURL": "{0}.bai".format(sorted_bam),
            "showSoftClips": "true",
        },
    ]
    with open(igv_json_file, "w") as handle:
        json.dump(track_config, handle)


def func(predict_tsv, fusion_db, oncokb_tsv, outd, samplename, output_bed_pe):
    # pydir = os.path.dirname(os.path.realpath(__file__))
    # dispose depth xls
    depth_df = pd.read_csv(
        f"{outd}/../{samplename}_depth.xls", sep="\t", index_col=None
    )
    uni_depth_dic = dict(zip(depth_df["Primer"], depth_df["Unique depth"]))
    # dispose coding effect tsv
    starfusion_coding_effect_df = pd.read_csv(
        f"{outd}/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
        sep="\t",
        index_col=None,
        # usecols=[0, 7, 9, 21],
    )
    coding_effect_dic = {
        f"{row['#FusionName']}_{row['LeftBreakpoint']}_{row['RightBreakpoint']}": row[
            "PROT_FUSION_TYPE"
        ]
        for _, row in starfusion_coding_effect_df.iterrows()
    }
    # dispose ctat introns file
    exon_skipping = read_ctat_splicing(
        f"{outd}/08_ctat.cancer.introns", f"{outd}/08_ctat.introns"
    )  # method
    # dispose transcript
    oncokb_df = pd.read_csv(oncokb_tsv, sep="\t", index_col=None, usecols=[5])
    oncokb_accession_lis = [
        i.split(".")[0] for i in oncokb_df["GRCh38 RefSeq"] if pd.notna(i)
    ]
    # match readid and primer name
    readid_primer_dic = match_readid_primer(outd)  # method

    # dispose predict tsv
    total_df = pd.DataFrame(
        columns=[
            "Fusion name",
            "Total count",
            "Alt count",
            "Fusion frequency%",
            "Fusion type",
            "LeftGene",
            "RightGene",
            "LeftBreakpoint",
            "RightBreakpoint",
            "LeftExonNumber",
            "RightExonNumber",
            "On target read count",
            "Off target read count",
            "Annotate",
        ]
    )
    all_fusion_intarget_read = []
    output_bedpe = open(output_bed_pe, "w")
    with open(predict_tsv, "r") as f:
        for line in f:
            newline = line.strip().split("\t")
            if re.match("#", line):
                continue
            fusion_name = newline[0]
            junction_read_count = int(newline[1])
            left_gene = newline[6].split("^")[0]
            left_break_point = newline[7]
            right_gene = newline[8].split("^")[0]
            right_break_point = newline[9]
            junction_readid = newline[10]
            annots = newline[18]
            transcript = fusion_name + "_" + left_break_point + "_" + right_break_point
            (
                primer_depth_dic,
                primer_readid_dic,
                fusion_intarget_read_lis,
                fusion_outoftarget_read_lis,
            ) = calculate_primer_frequency(junction_readid, readid_primer_dic)
            if not primer_depth_dic:
                continue
            fusion_all_primer_depth = 0
            for primer in primer_depth_dic:
                fusion_all_primer_depth += uni_depth_dic[primer]
            all_fusion_intarget_read += fusion_intarget_read_lis
            fusion_intarget_fre = round(
                len(fusion_intarget_read_lis) / fusion_all_primer_depth * 100, 2
            )
            left_primary_exon_number, left_exon_number = find_exon_number(
                left_break_point,
                f"{fusion_db}/ref_annot.gtf.sqldb",
                oncokb_accession_lis,
            )
            right_primary_exon_number, right_exon_number = find_exon_number(
                right_break_point,
                f"{fusion_db}/ref_annot.gtf.sqldb",
                oncokb_accession_lis,
            )
            if left_primary_exon_number != "NA" and right_primary_exon_number != "NA":
                fusion_name = "{0}({1})-{2}({3})".format(
                    left_gene,
                    left_primary_exon_number,
                    right_gene,
                    right_primary_exon_number,
                )
            else:
                fusion_name = "{0}-{1}".format(left_gene, right_gene)
            new_row = pd.DataFrame(
                {
                    "Fusion name": [fusion_name],
                    "Total count": [fusion_all_primer_depth],
                    "Alt count": [len(fusion_intarget_read_lis)],
                    "Fusion frequency%": [fusion_intarget_fre],
                    "Fusion type": [coding_effect_dic[transcript]],
                    "LeftGene": [left_gene],
                    "RightGene": [right_gene],
                    "LeftBreakpoint": [left_break_point],
                    "RightBreakpoint": [right_break_point],
                    "LeftExonNumber": [left_exon_number],
                    "RightExonNumber": [right_exon_number],
                    "On target read count": [len(fusion_intarget_read_lis)],
                    "Off target read count": [len(fusion_outoftarget_read_lis)],
                    "Annotate": [annots],
                }
            )
            total_df = pd.concat([total_df, new_row], ignore_index=True)
            # output bedpe
            if len(fusion_intarget_read_lis) < 5:
                continue
            if fusion_intarget_fre < 0.01:
                continue
            if fusion_all_primer_depth < 50:
                continue
            if coding_effect_dic[transcript] == "FRAMESHIFT":
                continue
            left_chr = left_break_point.split(":")[0]
            left_position = left_break_point.split(":")[1]
            right_chr = right_break_point.split(":")[0]
            right_position = right_break_point.split(":")[1]
            output_bedpe.write(
                "{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n".format(
                    left_chr, left_position, right_chr, right_position, fusion_name
                )
            )
    # check skipping fusion
    skipping = {}
    if exon_skipping:
        for i in exon_skipping:
            frequency = round(
                exon_skipping[i]["readnum"] / exon_skipping[i]["totalread"] * 100, 2
            )
            left_break_point = exon_skipping[i]["loc"].split("-")[0]
            right_break_point = (
                exon_skipping[i]["loc"].split(":")[0]
                + ":"
                + exon_skipping[i]["loc"].split("-")[-1]
            )
            left_chr = left_break_point.split(":")[0]
            left_position = left_break_point.split(":")[1]
            right_chr = right_break_point.split(":")[0]
            right_position = right_break_point.split(":")[1]
            if i == "MET_Exon14_Skipping":
                new_row = pd.DataFrame(
                    {
                        "Fusion name": ["MET_Exon14_Skipping"],
                        "Total count": [exon_skipping[i]["totalread"]],
                        "Alt count": [exon_skipping[i]["readnum"]],
                        "Fusion frequency%": [frequency],
                        "Fusion type": ["Non"],
                        "LeftGene": ["MET"],
                        "RightGene": ["MET"],
                        "LeftBreakpoint": [left_break_point],
                        "RightBreakpoint": [right_break_point],
                        "LeftExonNumber": ["Kit_Select:NM_000245:exon14"],
                        "RightExonNumber": ["Kit_Select:NM_000245:exon14"],
                        "On target read count": [exon_skipping[i]["readnum"]],
                        "Off target read count": [0],
                        "Annotate": ["Non"],
                    }
                )
                total_df = pd.concat([total_df, new_row], ignore_index=True)
                skipping["MET"] = {
                    "chr": left_chr,
                    "start": int(left_position),
                    "end": int(right_position),
                }
                output_bedpe.write(
                    "{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n".format(
                        left_chr, left_position, right_chr, right_position, i
                    )
                )
            elif i == "EGFR_Exon2-7_Skipping":
                new_row = pd.DataFrame(
                    {
                        "Fusion name": ["EGFR_Exon2-7_Skipping"],
                        "Total count": [exon_skipping[i]["totalread"]],
                        "Alt count": [exon_skipping[i]["readnum"]],
                        "Fusion frequency%": [frequency],
                        "Fusion type": ["Non"],
                        "LeftGene": ["EGFR"],
                        "RightGene": ["EGFR"],
                        "LeftBreakpoint": [left_break_point],
                        "RightBreakpoint": [right_break_point],
                        "LeftExonNumber": ["Kit_Select:NM_005228:exon2-7"],
                        "RightExonNumber": ["Kit_Select:NM_005228:exon2-7"],
                        "On target read count": [exon_skipping[i]["readnum"]],
                        "Off target read count": [0],
                        "Annotate": ["Non"],
                    }
                )
                total_df = pd.concat([total_df, new_row], ignore_index=True)
                skipping["EGFR"] = {
                    "chr": left_chr,
                    "start": int(left_position),
                    "end": int(right_position),
                }
                output_bedpe.write(
                    "{0}\t{1}\t{1}\t{2}\t{3}\t{3}\t{4}\n".format(
                        left_chr, left_position, right_chr, right_position, i
                    )
                )
    output_bedpe.close()
    extract_fusion_intarget_read_align(
        f"{outd}/07_star.Aligned.sortedByCoord.out.bam",
        f"{outd}/09_on_target_junction_reads.bam",
        all_fusion_intarget_read,
        f"{fusion_db}/ncbiRefSeq.sorted.txt",
        f"{outd}/09_trackConfigs.json",
        skipping,
    )  # method
    total_df = total_df.sort_values(by="Fusion frequency%", ascending=False)
    total_df.to_csv(f"{outd}/../{samplename}_fusion_total.xls", sep="\t", index=None)
    filter_df = total_df[
        (total_df["Total count"] >= 50)
        & (total_df["Alt count"] >= 5)
        & (total_df["Fusion frequency%"] >= 1)
        & (total_df["Fusion type"] != "FRAMESHIFT")
    ]
    filter_df.to_csv(f"{outd}/../{samplename}_fusion_filter.xls", sep="\t", index=None)

    # time_end = datetime.datetime.now()
    # print(
    #     f"\n[Time used]: \033[33m{round((time_end-time_start).total_seconds()/60, 2)}\033[0m min"
    # )


def main():
    argv = argparse_line()
    func(
        argv["predict_tsv"],
        argv["fusion_db"],
        argv["oncokb_tsv"],
        argv["outdir"],
        argv["prefix_name"],
        argv["output_bedpe"],
        # argv['threads'],
    )


if __name__ == "__main__":
    main()
