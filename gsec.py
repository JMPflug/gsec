#!/usr/bin/env python

import argparse
import os
import numpy as np
import pandas as pd
import collections
from functools import reduce

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


class Specimen:
    def __init__(self, in_file):
        self.in_file = in_file
        self.exclude_genes = []
        self.init_cov = None
        self.trim_cov = None

        self.sp_name = self.get_sp_name()

    def check_basecov(self, line_sample):
        self.line_sample = line_sample
        try:
            self.gene_id, self.pos, self.cov = self.line_sample.strip().split()
            if self.gene_id in self.exclude_genes:
                return False
            try:
                self.pos = int(self.pos)
            except ValueError:
                print("Second column in basecov file not an integer: {}".format(self.pos))
            try:
                self.cov = int(self.cov)
            except ValueError:
                print("Third column in basecov file not an integer: {}".format(self.cov))
        except ValueError:
            print("Check basecov file {}".format(self.in_file))
        return self.gene_id, self.pos, self.cov

    def get_sp_name(self):
        self.sp_name = os.path.split(self.in_file)[1]
        if "_basecov.txt" in self.sp_name:
            self.sp_name = self.sp_name.replace("_basecov.txt", "")
        return self.sp_name

    def initialize_loci(self):
        '''Opens basecov file and populates dict of genes and coverages,
        and renames, if necessary'''

        self.locus_dict = collections.OrderedDict()

        with open(self.in_file, 'r') as f:
            self.lines = [x for x in f if "#" not in x]
            for self.line in self.lines:
                self.gene_id, self.pos, self.cov = self.check_basecov(self.line)
                if self.gene_id:
                    self.populate_db(self.locus_dict, self.gene_id, self.pos, self.cov)

    def populate_db(self, locus_dict, gene_id, pos, cov):
        self.locus_dict = locus_dict
        self.id = id
        self.pos = pos
        self.cov = cov

        if self.gene_id not in self.locus_dict.keys():
            self.locus_dict[self.gene_id] = Locus(self.gene_id, self.pos, self.cov)
            # print("Init {} with pos {}, cov {}".format(self.gene_id, self.pos, self.cov))
        else:
            self.locus = self.locus_dict[self.gene_id]
            self.locus.pos.append(self.pos)
            self.locus.cov.append(self.cov)
            # print("adding to pos {}, cov {} to {}".format(self.pos, self.cov, self.gene_id))

    def trim_loci(self):
        for key, value in self.locus_dict.items():
            # print("Trimming {}".format(key))
            self.locus_dict[key].trim_bases(trim, mint)

    def sp_mean(self):
        self.all_covs = []
        self.col_name_init = "{}_Init".format(self.sp_name)
        self.col_name_trim = "{}".format(self.sp_name)

        for key, value in self.locus_dict.items():
            self.locus = self.locus_dict[key]
            # print(self.locus.init_mean, self.locus.trim_mean)
            self._dict = {'Locus': self.locus.id,
                          self.col_name_init: self.locus.init_mean,
                          self.col_name_trim: self.locus.trim_mean}
            # Remove initial coverage value
            del self._dict[self.col_name_init]
            self.all_covs.append(self._dict)

        return self.all_covs

class Locus:
    def __init__(self, id, pos, cov):
        self.id = id
        self.pos = [pos]
        self.cov = [cov]
        self.trim_cov = []
        self.init_mean = None
        self.trim_mean = None
        self.len_okay = True

    def set_means(self):
        """Calculates coverage means"""
        self.init_mean = np.mean(np.asarray(self.cov).astype(np.float))
        self.init_mean = round(self.init_mean, 4)
        self.trim_mean = np.mean(np.asarray(self.trim_cov).astype(np.float))
        self.trim_mean = round(self.trim_mean, 4)

    def trim_bases(self, trim, mint):
        self.trim = trim
        self.mint = mint
        self.cov_mod = self.cov[:]
        self.locus_len = len(self.pos)
        self.left_clip = []
        self.right_clip = []

        # Check length
        # print(self.locus_len, self.mint, self.trim)
        self.target_len = self.locus_len - self.mint - self.trim * 2
        if self.target_len <= 0:
            self.len_okay = False
            # print("Warning: Locus {} too short ({}). Minimum initial length is {}".format(self.id, str(self.locus_len), str(self.mint + self.trim * 2)))
        else:
            # Remove trim from each end of locus
            # print(self.cov)
            # Keeps trimmed portions. Currently unused.
            self.left_clip = []
            self.right_clip = []
            # Replace positions to be removed with *
            for i in range(0, self.trim):
                self.left_clip.append(self.cov[i])
                self.right_clip.append(self.cov[-(i+1)])
                self.cov_mod[i] = "*"
                self.cov_mod[-(i+1)] = "*"
            self.trim_cov = [y for y in self.cov_mod if y is not "*"]
            # print(self.trim_cov)
            self.set_means()
            # print(np.mean(np.asarray(self.trim_cov).astype(np.float)))

            # print("{} trimmed from {}bp to {}bp".format(self.id, self.locus_len, self.target_len))
            # print("\tCoverage went from {}bp to {}bp".format(self.init_mean, self.trim_mean))


def subset_by_iqr(df, column, iqr_mult=3):
    """Remove outliers from a dataframe by column, including optional
       whiskers, removing rows for which the column value are
       less than Q1-1.5IQR or greater than Q3+1.5IQR.
    Args:
        df (`:obj:pd.DataFrame`): A pandas dataframe to subset
        column (str): Name of the column to calculate the subset from.
        iqr_mult (float): Optional, loosen the IQR filter by a
                               factor of `iqr_mult` * IQR.
    Returns:
        (`:obj:pd.DataFrame`): Filtered dataframe
    """
    # Calculate Q1, Q2 and IQR
    q1 = df[column].quantile(0.25)
    q3 = df[column].quantile(0.75)
    iqr = q3 - q1
    # Apply filter with respect to IQR, including optional whiskers
    filter = (df[column] >= q1 - iqr_mult*iqr) & (df[column] <= q3 + iqr_mult*iqr)
    return df.loc[filter]


def iqr_trim(df, col):
    Q1 = df[col].quantile(0.25)
    Q3 = df[col].quantile(0.75)
    IQR = Q3 - Q1
    colst = str(col)
    filtered = df.query('(@Q1 - 1.5 * @IQR) <= @colst <= (@Q3 + 1.5 * @IQR)')
    # print(filtered)
    df.join(filtered, rsuffix='_filtered').boxplot()


def exclude_gene(self):
    '''Takes a list of gene names, separated by whitespace, and returns them as a list'''
    # If ignore_gene list is empty, return None
    if self.excluded_genes is None:
        self.exclude_genes = [None]
        return self.excluded_genes

    # Populate list of excluded genes
    self.excluded_genes = self.excluded_genes.split()
    return self.excluded_genes


def trim(covs_in, gene_id, trim_size, trim_fract):
    covs = covs_in[:]
    cov_len = len(covs_in)
    trim_covs = covs_in[trim_size:len(covs)-trim_size]
    #trim_covs_len = len(trim_covs)
    trim_covs_len = len(covs_in)
    trim_window = int(round(trim_covs_len*(1-trim_fract)))

    left_clip = []
    right_clip = []
    if trim_window/2 <= trim_size:
        for i in range(0, trim_size):
            left_clip.append(covs[i])
            right_clip.append(covs[-(i+1)])

            covs[i] = "*"
            covs[-(i+1)] = "*"
    else:
        for i in range(0, trim_window/2):
            left_clip.append(covs[i])
            right_clip.append(covs[-(i+1)])
            covs[i] = "*"
            covs[-(i+1)] = "*"

    window_covs = [y for y in covs if y is not "*"]
    win_mean = np.mean(window_covs)

    return window_covs


def main():
    parser = argparse.ArgumentParser(description="Removes a specified number of bases from the end of each gene in a BBMap base coverage file")

    parser.add_argument("-basecov", "-i",
                        help="Base coverage file.",
                        type=str)

    parser.add_argument("-exclude", "-e", default="trimmed_18S_MLSG_ trimmed_28S_MLSG_ trimmed_COI_MLSG_",
                        help="Optional list of gene names to ignore.",
                        type=str)

    parser.add_argument("-iqr", default=1.5,
                        help="IQR coefficient for excluding outliers.",
                        type=float)

    parser.add_argument("-trim", default=75,
                        help="Number of bases to remove from each end of locus.",
                        type=int)

    parser.add_argument("-min", default=100,
                        help="Minimum number of bases to retain a locus after trimming.",
                        type=int)

    parser.add_argument("-out", "-o",
                        help="Output name.", default="GSEC_out",
                        type=str)

    # parser.add_argument('--print-stats', dest='print_stats', action='store_true', help="Print detailed coverage information")

    # def check_fraction(value):
    #     fvalue = float(value)
    #     if fvalue > 1 or fvalue < 0:
    #         raise argparse.ArgumentTypeError("%s is an invalid fraction value. Enter a number between 0 and 1" % fvalue)
    #     return fvalue
    #
    # parser.add_argument("-trim_fraction", "-f", default=1, type=check_fraction,
    #                     help="Fraction of bases to retain from center of gene. Set\
    #                     to 1 to use all bases minus minimal terminal cutoff")

    # parser.set_defaults(print_stats=False)

    args = parser.parse_args()

    global trim_fraction, basecov, exclude, trim, iqr, mint

    trim_fraction = 1
    basecov = args.basecov
    exclude = args.exclude
    # print_stats = args.print_stats
    trim = args.trim
    iqr = args.iqr

    mint = args.min
    outname = args.out

    # trim_fraction = 1
    # basecov = "/Users/MaddisonLab/Documents/JMP/gsec_v0.3/gsec2/regier_baseco"
    # exclude = "trimmed_18S_MLSG_ trimmed_28S_MLSG_ trimmed_COI_MLSG_"
    # trim = 75
    # iqr = 1.5
    # mint = 100

    # in_basecovs = [x for x in os.listdir(basecov) if not (x.startswith(".")) and if os.path.isfile(x)]
    in_basecovs = [os.path.join(basecov, f) for f in os.listdir(basecov) if os.path.isfile(os.path.join(basecov, f)) and not f.startswith(".")]

    all_sp = []
    all_sp_orig = []

    for sp in in_basecovs:
        Sp = Specimen(sp)
        Sp.initialize_loci()
        Sp.trim_loci()
        df_or = pd.DataFrame.from_dict(Sp.sp_mean()).set_index('Locus')
        for _col in df_or.columns:
            df = subset_by_iqr(df_or, _col, iqr)

        all_sp.append(df)
        all_sp_orig.append(df_or)

    df_merged = reduce(lambda left,right: pd.merge(left,right,on=['Locus'],
                                                how='outer'), all_sp)

    df_means = df_merged.mean(skipna=True)

    print(df_means.to_csv(sep="\t", header=False).strip())

    outfile_name = outname + "_coverages.txt"
    with open(outfile_name, "w+") as f:
        f.write(df_means.to_csv(sep="\t", header=False))

    #Create boxplot

    df_or_merged = reduce(lambda left,right: pd.merge(left,right,on=['Locus'],
                                                how='outer'), all_sp_orig)
    df_graph = df_or_merged.melt(var_name='Specimen', value_name='Coverage')

    fig = px.box(df_graph, y='Coverage', x='Specimen', title="OrthoDB")

    fig.update_xaxes(title_font=dict(size=18, family='Helvetica'))
    fig.update_yaxes(title_font=dict(size=18, family='Helvetica'))
    fig.update_xaxes(tickangle=270, tickfont=dict(family='Helvetica', size=12))

    fig.update_layout(
        autosize=False,
        width=500,
        height=700,
        plot_bgcolor='rgba(0,0,0,0)',
        title=go.layout.Title(
        text="Locus Coverage",
        xref="paper",
        x=0.5
        )
        )

    fig.update_xaxes(automargin=True, linecolor="LightGray", mirror=True)
    fig.update_yaxes(automargin=True, gridcolor="LightGray", linecolor="LightGray", mirror=True)

    fig.update_traces(go.Box(marker=dict(
                                        outliercolor="Grey", color="Red", line=dict(outlierwidth=0.1, color="Green")),
                            line=dict(color="Red", width=1.2)),
                            line=dict(color="Gray")
                            )

    fig
    fig.write_image(outname + "_coverage_boxplot.pdf")


if __name__ == "__main__":
    main()














#
