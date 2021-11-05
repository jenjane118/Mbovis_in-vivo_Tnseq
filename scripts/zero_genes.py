#!/usr/bin python3
# coding=utf-8

""" Identify zero genes """

"""
Program:    zero_genes
File:       zero_genes.py
Version:    1.0
Date:       05.10.21
Function:   Create list of genes with zero insertions in both input/output files
Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/tn_seq

Institution:    Birkbeck University of London
                Project supervisors:  Dr. Irilenia Nobeli
                                  Dr. Sharon Kendall (RVC)
_____________________________________________________________________________
Description:
============
This program sums the reads at every insertion site in a gene over all the samples
and creates a list of genes that have nearly zero reads across all the insertion
sites in all the samples.

Usage:
======
zero_genes         SELF



"""

# ******************************************************************************
# Import libraries

import os, glob
import pandas as pd
import csv
import numpy as np

# ******************************************************************************

def make_wig_df(wig_file_dir, input_wig):
    """
    Function to open each wig file and add as column to dataframe of insertion sites.
    :param wig_file_dir: directory containing all wig files for analysis
    :return: df of insertion sites and number of reads in each sample
    (dim: columns are insertion sites (~73,000), rows are reads in each sample (27), first row is input lib)
    """

    with open(input_wig, 'r') as fi:
        input_table = pd.read_table(fi, sep=" ", skiprows=1, header=0, index_col=0)
    input_table.columns = ["reads"]
    # list of wig files
    dir = glob.glob(os.path.join(wig_file_dir, '*.wig'))
    # number of files in directory
    n_files = len(glob.glob(os.path.join(wig_file_dir, '*.wig')))
    wig_df = input_table
    for filename in dir:
        with open(filename, 'r') as f:
            sample_name = filename.replace("final_wigs_npremoved/perm_lung_tpp_", "")
        table = pd.read_table(filename, sep=" ", skiprows=1, header=0, index_col=0)
        table.columns = ["reads"]
        wig_df[sample_name] = table.reads

    return wig_df

# ******************************************************************************

def ta_list(tsv_file):
    """
    Function reads csv file of genes:insertion site coordinates and makes list.
    :param tsv_file: file of gene name and list of coordinates of insertions
    :return: list of genes and associated ta site coordinates
    """

    reader = csv.reader(open(tsv_file), delimiter = "\t")
    site_list = list(reader)
    return site_list


# ******************************************************************************

def zero_gene_list(wig_df, ta_list):
    """
    Function to iterate through wig dataframe, find sum of insertion sites indicated
    by site_list across all wig files (rows). Have one threshold for sum of reads at sites
    in one sample (5), and another threshold for sum of reads across samples at one site (55).

    :param wig_df: df that contains all the insertion sites and number of reads
                        for all the sample wig files
    :param site_list: list of gene and the list of coordinates for insertion sites
    :return: list of genes to be excluded from prot table
    """
    #test_genes = ["MB0004", "MB0047c", "MB0050", "MB0054", "MB0055", "MB0058", "MB0088", "MB0095"]
    #test_genes = ["MB1064c", "MB0521", "MB0471"]
    # make list of site counts and genes
    site_list = []  #list of site counts
    gene_list = []  #list of genes
    for gene, sites in ta_list[1:]:
        #if gene in test_genes:
        site_list.append(sites)
        gene_list.append(gene)
    zero_genes = []
    for i in range(0, len(gene_list)):
        sites = site_list[i].split(",")
        tester = site_tester(gene_list[i], sites, wig_df)
        if tester == True:
            zero_genes.append(gene_list[i])

    return zero_genes

# ******************************************************************************

def site_tester(gene, sites, wig_df):
    """
    Function to look at sum of reads at site across all samples and max reads in any one sample
    :param gene:
    :param sites:
    :param wig_df: dataframe of sites and reads in all samples
    :return: boolean operator T/F
    """

    test = True
    for coord in sites:
        coord = int(coord)         #convert counts to integers
        # check if sum of this site in all samples is below threshold
        site_sum = wig_df.loc[coord].sum()
        # maximum number of reads at site across samples
        max_reads = wig_df.loc[coord].max()
        if site_sum > 55 or max_reads > 5:
            test = False
            return test

    return test

########## main ###################################################################

if __name__ == "__main__":

    from datetime import datetime

    start = datetime.now()

    ta_file = "gene_ta_list.tsv"
    w_dir = "final_wigs_npremoved"
    i_wig = "perm_lung_tpp_MbA27.wig"

    l = ta_list(ta_file)
    
    t_arr = make_wig_df(w_dir, i_wig)

    zg = zero_gene_list(t_arr, l)
    print(len(zg))

    f = open("zero_genes.txt", "w")
    for gene in zg:
        f.write("%s\n" % gene)


    print(datetime.now() - start)
