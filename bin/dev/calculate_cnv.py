#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: zemin-chen@ebiotron.com

import numpy as np
import pandas as pd
import argparse, textwrap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D


def argparse_line():
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-cnr", metavar="FILE", help="input file, the cnvkit cnr file", required=True
    )
    parser.add_argument(
        "-cnv",
        metavar="FILE",
        help="output file, gene mean and median copy number",
        required=True,
    )
    parser.add_argument(
        "-png", metavar="PNG", help="output png, gene copy number plot", required=True
    )
    argv = vars(parser.parse_args())
    return argv


def clean_cnr_df(cnr_df):
    grouped = cnr_df.groupby("gene")
    cleaned_df = pd.DataFrame(columns=cnr_df.columns)
    for gene, df in grouped:
        q1 = df["copy_num"].quantile(0.25)
        q3 = df["copy_num"].quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        cleaned_df = pd.concat(
            [
                cleaned_df,
                df[(df["copy_num"] >= lower_bound) & (df["copy_num"] <= upper_bound)],
            ]
        )

    return cleaned_df


def calc_draw_cnv(cnr_file, cnv_file, png):
    cnr_df = pd.read_csv(cnr_file, header=0, sep="\t", index_col=False)

    # find CNV gene
    cnv_mask = cnr_df["gene"].str.contains(":")
    cnv_gene_series = cnr_df[cnv_mask]["gene"].str.split(":", expand=True)[0]

    cnv_gene_region_count = cnv_gene_series.value_counts()
    cnv_gene_region_count = cnv_gene_region_count[cnv_gene_region_count > 5]
    cnv_gene = cnv_gene_region_count.index.tolist()

    # remove CNV mark
    cnr_df["gene"] = cnr_df["gene"].str.split(":").str[0]

    # remove Antitarget region
    # cnr_df = cnr_df[cnr_df['gene'] != 'Antitarget']

    # calculate copy number
    cnr_df["copy_num"] = (2 ** (cnr_df["log2"])) * 2

    # remove gene outlier copy number
    cnr_df = clean_cnr_df(cnr_df)

    # calculate gene mean and median copy number
    gene_copy_num_df = (
        cnr_df.groupby("gene")["copy_num"].agg(["mean", "median"]).reset_index()
    )

    # calculate bin position
    cnr_df["position"] = cnr_df.apply(
        lambda row: int(row["start"] + (row["end"] - row["start"]) / 2), axis=1
    )

    # generate rank from position
    for chrom, df in cnr_df.groupby("chromosome"):
        df["position"] = df["position"].rank(method="dense")
        cnr_df.loc[df.index] = df

    # Create a Categorical data type with the desired sorting order
    chrom_order = pd.CategoricalDtype(
        [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
        ],
        ordered=True,
    )

    grouped = cnr_df.groupby("chromosome")
    sorted_grouped = sorted(grouped, key=lambda x: chrom_order.categories.get_loc(x[0]))

    # Create a dictionary to map each gene to a unique color
    np.random.seed(101)
    gene_colors = {}
    for gene in cnr_df.gene.unique():
        # Generate a random RGB color for each gene
        if gene in cnv_gene:
            gene_colors[gene] = np.random.rand(
                3,
            )
        else:
            gene_colors[gene] = (211 / 255, 211 / 255, 211 / 255, 0.6)

    chromosome_point = cnr_df.groupby("chromosome").size()
    chromosome_point = chromosome_point.sort_values(
        key=lambda x: x.index.astype(chrom_order)
    )
    # chromosome_point = np.log(chromosome_point)

    # Create grid layout
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(
        1, len(sorted_grouped), figure=fig, width_ratios=chromosome_point
    )

    max_copy_nums = int(
        cnr_df.loc[cnr_df["gene"].isin(cnv_gene)]
        .groupby("gene")["copy_num"]
        .max()
        .max()
    )

    for plot_num, (chr_num, group) in enumerate(sorted_grouped):
        ax = plt.subplot(gs[plot_num])
        ax.axhline(y=2, color="green", linestyle="--")
        ax.axhline(y=4, color="red", linestyle="--")
        for gene in group.gene.unique():
            color = gene_colors[gene]
            gene_data = group[group.gene == gene]
            ax.scatter(gene_data.position, gene_data.copy_num, color=color, s=15)
            if gene in cnv_gene:
                max_point = gene_data[
                    gene_data.copy_num == gene_data.copy_num.max()
                ].iloc[0]
                ax.annotate(
                    gene,
                    xy=(max_point.position, max_point.copy_num),
                    xytext=(-12, 10),
                    size=8,
                    textcoords="offset points",
                    va="center",
                )
                ###############################################################
                mean_copy_num = round(
                    gene_data.groupby("gene")["copy_num"].mean()[gene], 2
                )
                median_copy_num = round(
                    gene_data.groupby("gene")["copy_num"].median()[gene], 2
                )
                if mean_copy_num >= 3.5:
                    ax.axhline(y=mean_copy_num, color=gene_colors[gene], linestyle="--")
                    ax.axhline(
                        y=median_copy_num, color=gene_colors[gene], linestyle="--"
                    )
                    ax.annotate(
                        "mean="
                        + str(mean_copy_num)
                        + "\n"
                        + "median="
                        + str(median_copy_num),
                        xy=(max_point.position, max_point.copy_num),
                        xytext=(-50, 25),
                        size=8,
                        textcoords="offset points",
                        va="center",
                    )
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        margin = (group.position.max() - group.position.min()) / 10
        ax.set_xlim(group.position.min() - margin, group.position.max() + margin)
        if max_copy_nums < 10:
            plt.ylim(0, 10)
            plt.yticks(range(0, 11, 1))
        else:
            plt.ylim(0, max_copy_nums + 4)
            plt.yticks(range(0, max_copy_nums + 4, 1))
        ax.set_xlabel(chr_num)
        plt.setp(ax.get_xticklabels(), visible=False)
        if plot_num > 0:
            ax.tick_params(left=False)
            ax.spines["left"].set_visible(False)
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel("Copy number", fontsize=15)
        ax.autoscale_view()

    # Adjust subplot spacing
    plt.subplots_adjust(wspace=0)

    # Create a custom legend for the gene colors
    legend_elements = []
    for gene in cnv_gene:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=gene,
                markerfacecolor=gene_colors[gene],
                markersize=12,
            )
        )

    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 0.5), loc="center left")
    plt.suptitle(cnr_file.split(".")[0], fontsize=16, y=0.92)
    plt.savefig(png, dpi=300)

    gene_copy_num_df["mean"] = round(gene_copy_num_df["mean"], 3)
    gene_copy_num_df["median"] = round(gene_copy_num_df["median"], 3)
    gene_copy_num_df.rename(
        columns={
            "gene": "Gene",
            "mean": "CopyNumber(mean)",
            "median": "CopyNumber(median)",
        },
        inplace=True,
    )
    gene_copy_num_df.to_csv(cnv_file, sep="\t", index=False)


def main():
    argv = argparse_line()
    calc_draw_cnv(argv["cnr"], argv["cnv"], argv["png"])


if __name__ == "__main__":
    main()
