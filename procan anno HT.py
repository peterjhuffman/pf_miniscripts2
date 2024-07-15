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

CLx = ['Calu-6', 'ChaGo-K-1', 'HCC-78', 'NCI-H1648', 'NCI-H1650', 'SK-MES-1']


df = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan TAA vec_i.xlsx')
print('File read successfully.')



rnaq = True
protq = True
jimboq = True

for CL in CLx:
    if [x for x in cell_lines['pcells']].count(CL)>0:
        if protq:
            proteomics = pd.read_csv(path+f'proteomics\\{CL}_proteomics.csv')
            proteomics = proteomics[proteomics.columns[2:]].T.reset_index()
            del proteomics[0]
            proteomics.columns = ['Accession', f'ProCAN Abundance: {CL}']
            proteomics=proteomics.fillna(0)
            df = pd.merge(df,proteomics, how='left', on='Accession')
            df[f'ProCAN Abundance: {CL}'] = df[f'ProCAN Abundance: {CL}'].fillna(0)
            df[f'ProCAN Abundance: {CL}'] = df[f'ProCAN Abundance: {CL}'].astype(float)
    else:
        print('Cell line not aquired in Proteomics.')
    
    
    if [x for x in cell_lines['rcells']].count(CL)>0:
        if rnaq:
            rnaseq = pd.read_csv(path+f'rnaseq tpm\\{CL}_rnatpm.csv')
            rnaseq.columns = [f'ProCAN TPM: {CL}', 'Gene Symbol']
            df = pd.merge(df, rnaseq, how='left',on='Gene Symbol')
        if [x for x in cell_lines['pcells']].count(CL)>0:
            if jimboq:
                df['ProCAN Score d1'] = 0
                df['ProCAN Score d2'] = 0
                tpm_mask = df[f'ProCAN TPM: {CL}']>1
                prot_mask = df[f'ProCAN Abundance: {CL}']>0
                df.loc[tpm_mask, 'ProCAN Score d1'] = 1
                df.loc[prot_mask, 'ProCAN Score d2'] = 1
                df['ProCAN AcqScore'] = df['ProCAN Score d1']+df['ProCAN Score d2']
                del df['ProCAN Score d1']
                del df['ProCAN Score d2']
        else:
            print('Cell line not aquired in Proteomics.')
    else:
        print('Cell line not aquired in RNAseq.')




df.to_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan TAA vec_proCanno.xlsx',index=False)
print('File saved')


input('Procan ANNO complete. press any key to exit.')














