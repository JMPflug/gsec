#!/usr/bin/env python

import argparse
import os
import numpy as np
import collections

parser = argparse.ArgumentParser(description="Removes a specified number of bases from the end of each gene in a BBMap base coverage file")

parser.add_argument("-basecov", "-i",
                    help="Base coverage file.",
                    type=str)

parser.add_argument("-min_trim_size", "-s",
                    help="Minimum number of bases to trim from both ends.",
                    type=int)

parser.add_argument("-trans_table", "-t", default=None,
                    help="Optional table of names to convert.",
                    type=str)

parser.add_argument("-basecov_text", "-b", default="basecov",
                    help="Text string used to indicate basecov files if a directory is given as the input for -basecov.",
                    type=str)

parser.add_argument("-exclude_gene", "-e", default="trimmed_18S_MLSG_ trimmed_28S_MLSG_ trimmed_COI_MLSG_",
                    help="Optional list of gene names to ignore.",
                    type=str)

parser.add_argument("-read_length", "-l", default=None,
                    help="",
                    type=int)

parser.add_argument("-read_number", "-r", default=None,
                    help="",
                    type=int)

parser.add_argument("-read_stats", "-d", default=None,
                    help="A tab delimited file of library names, number of reads, and read lengths. Optional for calculating GSE directly",
                    type=str)

parser.add_argument("-iqr", default=1.5,
                    help="IQR coefficient for excluding outliers. If --exclude-outliers flag not set, this is ignored",
                    type=float)

parser.add_argument('--print-cov', dest='print_cov', action='store_true', help="Print detailed coverage information")

parser.add_argument('--print-basecov', dest='print_basecov', action='store_true', help="Print long form base coverage")

parser.add_argument('--print-stats', dest='print_stats', action='store_true', help="Suppress stats printout")

parser.add_argument('--inexact-match', dest='inexact_match', action='store_true', help="Allow exclude genes to support inexact matches. Useful for excluding groups of genes containing particular substrings.")

parser.add_argument('--exclude-outliers', dest='exclude_outliers', action='store_true', help="Excludes loci based on IQR. Adjust IQR coefficient usint -iqr parameter")

parser.add_argument('--exclude-long', dest='exclude_long', action='store_true', help="Excludes loci that are unusually long, which may cause inflated PerBase estimates")

parser.set_defaults(print_stats=False, print_cov=False, print_basecov=False, inexact_match=False, exclude_outliers=False, exclude_long=False)


def check_fraction(value):
    fvalue = float(value)
    if fvalue > 1 or fvalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid fraction value. Enter a number between 0 and 1" % fvalue)
    return fvalue


parser.add_argument("-trim_fraction", "-f", default=1, type=check_fraction,
                    help="Fraction of bases to retain from center of gene. Set\
                    to 1 to use all bases minus minimal terminal cutoff")

args = parser.parse_args()

in_basecov = args.basecov
trim_size = args.min_trim_size
trans_table = args.trans_table
ignore_gene = args.exclude_gene
trim_fraction = args.trim_fraction
read_number = args.read_number
read_length = args.read_length
read_stats = args.read_stats
iqr_coef = args.iqr

inexact_match = args.inexact_match
print_trim_covs = args.print_cov
print_basecov = args.print_basecov
print_stats = args.print_stats
basecov_text = args.basecov_text
exclude_outliers = args.exclude_outliers
length_outliers = args.exclude_long

perbase = True

#./trim_coverage_v0.9.py -s 75 -i /Volumes/Genomics_Main/UCE_Beetle_Collaboration_Assemblies/Reads/filtered_reads/Regier_Mapping_Chlaenius_sericeus_v2/Chlaenius_sericeus_DNA4821_basecov.txt  -t /Volumes/Genomics_Main/UCE_Beetle_Collaboration_Assemblies/Reads/filtered_reads/scripts/regier_id_trans_table_v3.txt --inexact-match -d read_stats.txt

def trans_names():
    '''Takes a list containing a pair of old and new gene names, one per line,
    and renames the corresponding old name to the new one'''

    # If no list is given, set tnames_dict to None and return
    if trans_table is None:
        tnames_dict = None
        return tnames_dict
    tnames_dict = {}
    # Open table and populate tname_dict
    with open(trans_table, 'r') as g:
        for line in g:
            old_name, new_name = line.split()
            tnames_dict[old_name] = new_name
    return tnames_dict


def exclude_gene():
    '''Takes a list of gene names, separated by whitespace, and returns them as a list'''
    # If ignore_gene list is empty, return None
    if ignore_gene is None:
        exclude_genes = [None]
        return exclude_genes
    ignore_genes = ignore_gene.split()
    exclude_genes = []
    # Populate list of excluded genes
    if ignore_genes is not None:
        for gene in ignore_genes:
            exclude_genes.append(gene)
    return exclude_genes


def populate_and_rename(in_file):
    '''Opens basecov file and populates dict of genes and coverages,
    and renames, if necessary'''

    tnames_dict = trans_names()
    # for key in tnames_dict:
    #     print key
    exclude_genes = exclude_gene()
    pos_dict = collections.OrderedDict()

    with open(in_file, 'r') as f:
        for line in f:
            include = True
            # Ignore header line
            if "#" not in line:
                gene_id, pos, cov = line.strip().split()
                # Rename gene_id if found in tname_dict
                if tnames_dict is not None and gene_id in tnames_dict:
                    gene_id = tnames_dict[gene_id]
                cov = int(cov)
                if gene_id not in pos_dict and gene_id not in exclude_genes:
                    for exc_gene in exclude_genes:
                        if exc_gene.strip() in gene_id.strip():
                            if inexact_match is True:
                                include = False
                    if include is True:
                        pos_dict[gene_id] = [cov]
                        #print gene_id

                elif gene_id in pos_dict and gene_id not in exclude_genes and include is True:
                    for exc_gene in exclude_genes:
                        if exc_gene.strip() in gene_id.strip():
                            print(exc_gene, gene_id)
                        else:
                            pos_dict[gene_id].append(cov)

    return pos_dict


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

    if print_basecov is True: 
        pos = 0
        for cov in covs_in:
            print("{}\t{}\t{}".format(gene_id, pos, cov))
            pos = pos+1

        for tcov in covs:
            print("{}\t{}\t{}".format(gene_id+"_trimmed", pos, tcov))
            pos = pos+1


    window_covs = [y for y in covs if y is not "*"]
    win_mean = np.mean(window_covs)

    if print_trim_covs is True:
        left_mean = np.mean(left_clip)
        right_mean = np.mean(right_clip)
        print(gene_id)
        print("\tTrimmed {}\tLen: {} Center: {} Fraction: {} Trim Window: {} Rep Count: {}".format(gene_id, cov_len, cov_len/2, trim_fract, trim_window, covs.count("*"))) 
        print("\tMeans Window: {} Left: {} Right: {}".format(win_mean, left_mean, right_mean))

    return window_covs


def print_pos():
    entry_list = []
    pos_dict = populate_and_rename()
    for key in pos_dict:
        covs = pos_dict[key]
        #print covs
        #trim_covs = covs[trim_size-1:len(covs)-trim_size-1]
        if len(covs) > trim_size*2:
            #trim_covs = covs[trim_size:len(covs)-trim_size]
            trim_covs = trim(covs, key)
            # print "trim_covs"
            # print trim_covs
        else:
            trim_covs = [0]
        #entry = '{}\t{}\t{}\t{}\t{}'.format(key, len(covs), len(trim_covs), np.mean(covs), np.mean(trim_covs))
        lenxcov = len(covs)*np.mean(covs)
        trimlenxcov = len(trim_covs)*np.mean(trim_covs)
        entry = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(key, len(covs), len(trim_covs), np.mean(covs), np.mean(trim_covs), lenxcov, trimlenxcov)
        if key is not None:
            entry_list.append(entry)

    seq_name = in_basecov.split('/')[-1]
    
    if print_stats is True:
        #print "\n"
        print('\t\t\t{}'.format(seq_name))
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format("#ID", "Length", "TrimLen", "Coverage", "TrimCov", "Len*Cov", "Trim Len*Cov"))

        for entry in entry_list:
            print(entry)

def print_pos_varpar(in_file, trim_sizes, trim_fracts):
    # trim_size = args.min_trim_size
    # trim_fraction = args.trim_fraction

    entry_list = []
    pos_dict = populate_and_rename(in_file)
    len_cut, len_median = find_length_outliers(pos_dict)
    len_warn = False
    len_warn_list = []
    header_list = ["#ID", "Length_Init", "Cov_Init", "PerBase_Init"]
    for temp_trim_size in trim_sizes:
        for temp_trim_fract in trim_fracts:
            len_name = "Length_Tr{}Fr{}\tCov_Tr{}Fr{}\tPerBase_Tr{}Fr{}".format(temp_trim_size, temp_trim_fract, temp_trim_size, temp_trim_fract, temp_trim_size, temp_trim_fract)
            len_name = len_name.split("\t")
            header_list.extend(len_name)
    for key in pos_dict:
        covs = pos_dict[key]
        raw_mean_cov = np.mean(covs)
        init_entry = [key, str(len(covs)), str(raw_mean_cov), str(len(covs)*raw_mean_cov)]
        if len(covs) > len_cut and length_outliers is True:
            print("WARNING: {} in {} is unusually large ({} bases; the median is {} bases)".format(key, in_file, len(covs), len_median))
            continue
        elif len(covs) > len_cut and length_outliers is False:
            len_warn = True
            len_warn_list.append([key, len(covs)])
        entry_list.append(init_entry)
        trim_sizes_entries = []
        #trim_sizes_entries.append(in_file.split('/')[-1])
        for temp_trim_size in trim_sizes:
            for temp_trim_fract in trim_fracts:
                trim_size = temp_trim_size
                #len_name = "Lenght Trim{}".format(trim_size)
                #header_list.append(len_name)
                if len(covs) > trim_size*2:
                    trim_cov = trim(covs, key, trim_size, temp_trim_fract)
                else:
                    trim_cov = [0]
                #print key, str(np.mean(trim_cov))
                trim_sizes_entries.append( str(len(trim_cov)) )
                trim_sizes_entries.append( str(np.mean(trim_cov)) )
                trim_sizes_entries.append( str(np.mean(trim_cov)*len(trim_cov) ))

        init_entry.extend(trim_sizes_entries)

    seq_name = in_file.split('/')[-1]

    per_file_entry = []
    per_file_entry.append(header_list)

    #print '\t\t\t{}'.format(seq_name)

    if len_warn is True:
        print("WARNING: Some loci are unusually long, which may result in inflated PerBase coverages. Consider using the \'--exclude-long\' flag\nAbnormal loci:")
        for l in len_warn_list:
            print("\t{} ({} bases)".format(l[0], l[1]))
    if print_stats is True:
        print("\n{} Mapping Statistics".format(seq_name))
        print('\t'.join(header_list))

    for entry in entry_list:
        if print_stats is True:
            print('\t'.join(entry))
        per_file_entry.append(entry)

    per_file_array = np.array(per_file_entry)

    return per_file_array


def find_length_outliers(pos_dict):
    len_list = []
    for key in pos_dict:
        len_list.append(len(pos_dict[key]))
    len_list = np.array(len_list)
    iqr_screen(len_list, "")
    q1, q3 = np.percentile(len_list, [25,75])
    len_cut = q3
    return len_cut, np.median(len_list)


def print_summary(in_array):
    print("{} Summary".format(in_basecov.split("/")[-1]))

    if read_length is not None and read_number is not None:
        print("Genome Size Estimate")
        print("Field\tGSE")
        for x in in_array.T[1:]:
            col = x[1:]
            col = col.astype('f')
            if "Cov" in x[0]:
                print("{}\t{}".format(x[0], (read_length*read_number)/np.mean(col) ))
    else:
        print("Average Coverage")
        print("Field\tMean")
        for x in in_array.T[1:]:
            col = x[1:]
            col = col.astype('f')
            if "Cov" in x[0]:
                print("{}\t{}".format(x[0], np.mean(col)))


def dir_print_summary(in_array, seq_name):
    ret_list = []
    read_stat_list = []
    if read_stats is not None:
        read_stats_dict = parse_read_stats()
        #print "GSE\t{}".format(seq_name)
        ret_list.append(["GSE", seq_name])
        for x in in_array.T[1:]:
            col = x[1:]
            col = col.astype('f')
            if "Cov" in x[0]:

                for key in read_stats_dict:
                    if key in seq_name:
                        read_stat_list = read_stats_dict[key]
                        read_length = int(read_stat_list[0])
                        read_number = int(read_stat_list[1])
                if not read_stat_list:
                    raise ValueError('No statistics found for ' + seq_name)


                gse = (read_number*read_length)/np.mean(col)

                #print "{}\t{}".format(x[0], gse)
                ret_list.append([x[0], str(gse)])

        if perbase is True:
            for i, x in enumerate(in_array.T[1:]):
                col = x[1:]
                col = col.astype('f')
                if "PerBase" in x[0]:

                    for key in read_stats_dict:
                        if key in seq_name:
                            read_stat_list = read_stats_dict[key]
                            read_length = int(read_stat_list[0])
                            read_number = int(read_stat_list[1])
                    if not read_stat_list:
                        raise ValueError('No statistics found for ' + seq_name)

                    len_col = in_array[1:,i-1]
                    len_col_int = [ int(y) for y in len_col ]
                    #print x[0], sum(col), in_array[0,i-1], sum(len_col_int), (sum(col)/sum(len_col_int))
                    perbase_cov = sum(col)/sum(len_col_int)

                    gse = (read_number*read_length)/perbase_cov

                    ret_list.append( [x[0], str(gse)])

    else:
        ret_list.append(["Coverage", seq_name])
        for x in in_array.T[1:]:
            col = x[1:]
            col = col.astype('f')

            if "Cov" in x[0]:
                col = sorted(col)
                col_per = np.array(col)
                col_per = np.asfarray(col_per, float)
                q1, q3 = np.percentile(col_per, [25,75])
                #print x[0], q1, q3
                iqr = q3 - q1
                #print x[0]
                #print "Q3", q3, "Q1", q1, "IQR", iqr
                lower_bound = q1 - (iqr_coef * iqr)
                upper_bound = q3 + (iqr_coef * iqr)
                #print "\n"+str(col_per)+"\n"
                count = 0
                for n in col:
                    if exclude_outliers is True and "_Init" not in x[0] and (n < lower_bound or n > upper_bound):
                        #print col[count]
                        del col[count]
                        #print "Outlier found at IQR={}: {} {} {}\n{} - {}".format(iqr_coef, seq_name, x[0], n, lower_bound, upper_bound)
                    count = count + 1
                #print "{}\t{}\t{}".format(x[0], np.mean(col), np.median(col))
                ret_list.append([x[0], str(np.mean(col).item())])
    
        if perbase is True:
            for i, x in enumerate(in_array.T[1:]):
                col = x[1:]
                col = col.astype('f')
                if "PerBase" in x[0]:
                    if exclude_outliers is True and "_Init" not in x[0]:
                        col = iqr_screen(col, seq_name)
                    len_col = in_array[1:,i-1]
                    len_col_int = [int(y) for y in len_col]
                    #print x[0], sum(col), in_array[0,i-1], sum(len_col_int), (sum(col)/sum(len_col_int))

                    perbase_cov = sum(col)/sum(len_col_int)
                    ret_list.append( [x[0], str(perbase_cov)])


    return ret_list

def iqr_screen(col, seq_name):
    col_per = np.array(col)
    col_per = np.asfarray(col_per, float)
    q1, q3 = np.percentile(col_per, [25,75])
    iqr = q3 - q1

    #print "Q3", q3, "Q1", q1, "IQR", iqr
    lower_bound = q1 - (iqr_coef * iqr)
    upper_bound = q3 + (iqr_coef * iqr)
    #print "\n"+str(col_per)+"\n"
    count = 0
    col = np.ndarray.tolist(col)
    for n in col:
        if exclude_outliers is True and (n < lower_bound or n > upper_bound):
            #print col[count]
            del col[count]
            #print "Outlier found at IQR={}: {} {}\n{} - {}".format(iqr_coef, seq_name, n, lower_bound, upper_bound)
        count = count + 1
    
    return np.array(col)

def parse_read_stats():
    read_stats_dict = {}
    with open(read_stats, "r") as f:
        for line in f:
            if "#" not in line:
                lib_name, read_number, read_length = line.split()
                read_stats_dict[lib_name] = [read_length, read_number]
    return read_stats_dict


def slistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]


def main():
    if os.path.isdir(in_basecov) is True and os.path.isfile(in_basecov) is False:
        print("Processing directory: {}".format(in_basecov))
        cov_summaries = []
        for in_file in slistdir(in_basecov):
            if basecov_text in in_file:
                in_file = os.path.join(in_basecov, in_file)
                seq_name = in_file.split('/')[-1].split("_" + basecov_text)[0]
                output = print_pos_varpar(in_file, [75], [.5])
                #output = print_pos_varpar(in_file, [75], [1, .5])

                temp_covs = dir_print_summary(output, seq_name)
                if any(cov_summaries) is False:
                    for cov in temp_covs:
                        cov_summaries.append(cov)
                else:
                    for i, cov in enumerate(temp_covs):
                        cov_summaries[i].append(str(cov[1]))
        # print "\n"
        # print in_basecov.split('/')[-2]
        # for line in cov_summaries:
        #     print "\t".join(line)

        print("\n")
        tp_cov_sum = list(map(list, list(zip(*cov_summaries))))
        for line in tp_cov_sum:
            print("\t".join(line))

    
    else:
        for in_file in in_basecov.split():
            print("Processing file: {}".format(in_basecov))
            seq_name = in_file
            #output = print_pos_varpar(in_file, [75, 90], [1, .9, .8, .7, .6, .5, .4, .3, .2, .1])
            output = print_pos_varpar(in_file, [90], [1, .5])

            print_summary(output)

    print("\nDone  ")




main()













