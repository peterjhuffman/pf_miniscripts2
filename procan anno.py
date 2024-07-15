# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:06:24 2024

@author: HUFFMP
"""

import pandas as pd
import numpy as np
import sys
from time import sleep
import pathlib

def crash():
    print('\nExiting procan anno. \n\n')
    sleep(5.2)
    sys.exit()


path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\procan slices\\'

# print('Loading ProCAN databases...')
# proteomics = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/general - chemoproteomics/procan proteomic scores.xlsx')
# print('proteomics read.')
# rnaseq = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/general - chemoproteomics/procan_rnaseq/procan rnaseq tpm.xlsx')
# print('rnaseq tpm read.')

cell_lines = pd.read_excel(path+'cell_lines.xlsx')

switch = True
while switch:
    CL = str(input('Cell Line? : ').upper())

    switch = ([x for x in cell_lines['pcells']].count(CL)+[x for x in cell_lines['rcells']].count(CL))==0
    

    # switch = False
    if switch:
        print('Cell Line not recognized.')


if [x for x in cell_lines['pcells']].count(CL)>0:
    print('Cell Line acquired in Proteomics.')
if [x for x in cell_lines['rcells']].count(CL)>0:
    print('Cell Line acquired in Transcriptomics.')


switch2 = True
while switch2:
    filename = input('Protein Export Filename? : ')
    if filename[0]=='f':
        filename = filename[5:]
    fileSuffix = pathlib.Path(filename).suffix
    foldersav = pathlib.Path(filename).stem
    filetrack = filename[:-(len(foldersav)+len(fileSuffix))]
    print('Searching for and reading file.')
    if fileSuffix == '.xlsx':
        df = pd.read_excel(filename)
        switch2=False
    elif fileSuffix == '.csv':
        df = pd.read_csv(filename)
        switch2=False
    else:
        print('File type not recognized.')
    if (not switch2):
        if [x for x in df.columns].count('Accession')==0:
            switch2=True
            print('Your file is missing accession. add that its important')
        if [x for x in df.columns].count('Gene Symbol')==0:
            switch2=True
            print('Your file is missing gene symbol. add that its important')

print('File read successfully.')



opps = ['Y', 'N']

rnaq = str(input('Add RNAseq data? (Y/N):  ')).upper()
while opps.count(rnaq)==0:
    print('Type Y for yes or N for no.')
    rnaq = str(input('Add RNAseq data? (Y/N):  ')).upper()

protq = str(input('Add Proteomics data? (Y/N):  ')).upper()
while opps.count(protq)==0:
    print('Type Y for yes or N for no.')
    protq = str(input('Add Proteomics data? (Y/N):  ')).upper()

jimboq = str(input('Add multiomic summary? (Y/N):  ')).upper()
while opps.count(protq)==0:
    print('Type Y for yes or N for no.')
    jimboq = str(input('Add multiomic summary? (Y/N):  ')).upper()



rnaq = rnaq=='Y'
protq = protq=='Y'
jimboq = jimboq=='Y'

if [x for x in cell_lines['pcells']].count(CL)>0:
    if protq:
        proteomics = pd.read_csv(path+f'proteomics\\{CL}_proteomics.csv')
        proteomics = proteomics[proteomics.columns[2:]].T.reset_index()
        del proteomics[0]
        proteomics.columns = ['Accession', 'ProCAN Proteomic Abundance']
        proteomics=proteomics.fillna(0)
        df = pd.merge(df,proteomics, how='left', on='Accession')
        df['ProCAN Proteomic Abundance'] = df['ProCAN Proteomic Abundance'].fillna(0)
        df['ProCAN Proteomic Abundance'] = df['ProCAN Proteomic Abundance'].astype(float)
else:
    print('Cell line not aquired in Proteomics.')


if [x for x in cell_lines['rcells']].count(CL)>0:
    if rnaq:
        rnaseq = pd.read_csv(path+f'rnaseq tpm\\{CL}_rnatpm.csv')
        rnaseq.columns = ['ProCAN RNAseq TPM', 'Gene Symbol']
        df = pd.merge(df, rnaseq, how='left',on='Gene Symbol')
    if [x for x in cell_lines['pcells']].count(CL)>0:
        if jimboq:
            df['ProCAN Score d1'] = 0
            df['ProCAN Score d2'] = 0
            tpm_mask = df['ProCAN RNAseq TPM']>1
            prot_mask = df['ProCAN Proteomic Abundance']>0
            df.loc[tpm_mask, 'ProCAN Score d1'] = 1
            df.loc[prot_mask, 'ProCAN Score d2'] = 1
            df['ProCAN AcqScore'] = df['ProCAN Score d1']+df['ProCAN Score d2']
            del df['ProCAN Score d1']
            del df['ProCAN Score d2']
    else:
        print('Cell line not aquired in Proteomics.')
else:
    print('Cell line not aquired in RNAseq.')




df['ProCAN Proteomic Abundance']=df['ProCAN Proteomic Abundance'].astype(float)
df.to_excel(filetrack+foldersav+'_proCanno.xlsx',index=False)
print('File saved')


input('Procan ANNO complete. press any key to exit.')















