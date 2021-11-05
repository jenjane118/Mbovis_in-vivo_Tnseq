#!/usr/bin python3
# coding=utf-8

""" Remove non-permissive sites """

"""
Program:    non_permissive
File:       non_permissive.py
Version:    1.0
Date:       24.04.21
Function:   Finds location of non-permissive sequence motifs in specified genome.
Eliminates non-permissive TA sites in genome.

Author:     Jennifer J. Stiens
            j.j.stiens@gmail.com
            https://github.com/jenjane118/tn_seq

Institution:    Birkbeck University of London
                Project supervisors:  Dr. Irilenia Nobeli
                                  Dr. Sharon Kendall (RVC)
_____________________________________________________________________________
Description:
============
This program uses the sequence motif found to be non-permissive for Himar1 insertions,
(as described by DeJesus et al, 2017) to find non-permissive TA sites in specified
genome and remove the identified sites from the insertion files (.wig).

Usage:
======
non_permissive         SELF


"""

# **********************************************************************************

def open_ta_list(file):
    """
    Read the .csv file of non-permissive sites and make into list of non-permissive TA sites.

    Input                   file            file of non-perm sites
    Output                  np_sites        list of non-perm sites

    """

    import pandas as pd

    df = pd.read_csv(file, header=1, low_memory=False)
    df = df[['Coordinate', 'Matches Non Permissive Motif']]
    np_sites = []
    for i in range(len(df)):
        if df['Matches Non Permissive Motif'][i] == True:
            np_sites.append(df['Coordinate'][i])

    return (np_sites)


###################################################################################

def remove_np_sites(wig_file, np_sites):
    """
    Open .wig file of insertion sites, delete non-permissive sites, save new .wig file

    Input               wig_file            file of insertion sites (two columns, position/no insertions)
                        np_sites            python list of positions of non-permissive sites
    Output              new_wig_file        new file of insertion sites (all permissive)

    """
    import pandas as pd
    import os

    # open file and record comment line
    f = open(os.path.expanduser(wig_file), 'r')
    comment = f.readline()
    df = pd.read_csv(wig_file, comment="#", sep=" ", header=0)
    f.close()
    # create name of new file
    new_filename = wig_file.split("/")[-1]
    new_wig_file = "~/git/tn_seq/lung_tnseq/tpp_genewiz1/np_removed/perm_" + new_filename
    # open outfile and put in comment
    outfile = open(os.path.expanduser(new_wig_file), 'w')
    outfile.write(comment)
    # find sites not in non-permissable list and write to outfile
    perm_sites = df[~df['variableStep'].isin(np_sites)]
    perm_sites.to_csv(outfile, index=False, sep=' ')
    outfile.close()

    return outfile


###################################################################################

def find_np_sites(motif, offset, fasta):
    """ Use motif finder to find all non-permissive TA sites in specified genome
    Input       motif               motif for non-permissive sites
                offset              number of nt before TA in motif
                fasta               fasta file for specified genome
    Output      np_sites            list of non-permissive positions

    """
    import os
    from Bio import SeqIO
    from Bio import SeqUtils

    # parse fastq file to extract sequence and convert to string
    fasta_file = os.path.expanduser(fasta)
    for record in SeqIO.parse(fasta_file, "fasta"):
        search_seq = str(record.seq)
    #search_seq = fasta
    # find matches for non-permissive motif
    matches = SeqUtils.nt_search(search_seq.upper(), motif)
    # add one to each position to make index=1 to match positions in wig file
    indexed_positions = []
    for i in range(1, len(matches)):
        new_position = matches[i] + 1 + offset
        indexed_positions.append(new_position)
    # create name of new file
    np_file = "~/git/tn_seq/lung_tnseq/bovis_np_sites.txt"
    # open outfile and put in comment
    outfile = open(os.path.expanduser(np_file), 'w')
    outfile.write('\n'.join(map(str, indexed_positions)))
    outfile.close()
    return indexed_positions


###################################################################################

def test_tb():
    # use published list of non-permissive motif locations
    no_sites = open_ta_list('~/tn_seq/data/dejesus/supp_2_non-permissive.csv')

    print(len(no_sites))

    a = remove_np_sites('~/Data/Dejesus_wig/SRR411328_1.wig', no_sites)

    print(a)


# **********************************************************************************


def test_fasta():

    np_motif_med = 'SGNTANCS'
    np_motif_long = 'GCGNTANCGC'
    bovis_fasta = "Mbovis_AF2122_97.fasta"
    tb_fasta = "~/git/mtb_modules/ref_seqs/Mtb_H37Rv.fasta"
    # create list of non-permissive positions
    no_sites = find_np_sites(np_motif_med, 3, bovis_fasta)
    # remove sites from insertion file, write new file
    new_file = remove_np_sites('~/git/tn_seq/lung_tnseq/tpp_genewiz1/lung_tpp_MbA013.wig', no_sites)
    print(new_file)


########## main ###################################################################

if __name__ == "__main__":

    test_fasta()
