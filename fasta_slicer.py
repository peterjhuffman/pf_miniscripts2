# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:54:20 2023

@author: HUFFMP
"""

"""

The purpose of this script is to take two files and return a truncated version of file 2.

File 1 is a dataframe that contains a list of accession IDs.
File 2 is a fasta file of the whole human proteome.

This program returns a fasta file containing only the proteins with accession IDs in file 1.

"""

import pandas as pd
import numpy as np

matlas = pd.read_csv('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\membrane proteomics\\'+
                     'NatCancer2021_surfaceProteome_tableS2_illustrated.csv')

matCodes = [str(acc) for acc in matlas['accession'].unique()]
fasta = open("\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\homo sapiens.fasta", "r")

fastalines = fasta.readlines()
new_fasta = []
switch = False

for line in fastalines:
    if line[0]== '>':
        res = [acc for acc in matCodes if(acc in line)]
        if bool(res):
            switch = True
        else:
            switch = False
    if switch:
        new_fasta.append(line)

fasta.close()

f = open("\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\homosapiens_natcancer_membraneatlas_rts.fasta", "a")
f.writelines(new_fasta)
f.close()

