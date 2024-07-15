# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:34:40 2024

@author: HUFFMP
"""
import pandas as pd
import numpy as np

filepath ='\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\ubr5global\\'
file = 'PFEC_PH_20240301_ubr5glob_ms_analysis.xlsx'

file1 = open(filepath+'prediction.txt', 'r')
xprediction = file1.readlines()[1:]
file1.close()

file2 = open(filepath+'known.txt', 'r')
xknown = file2.readlines()[1:]
file2.close()

known = [x.split('\t')[1] for x in xknown]
prediction = [x.split('\t')[3] for x in xprediction]

df = pd.read_excel(filepath+file)

# xf.loc[bhc_mask&s0_mask, f'{cat_code} Label FDR'] = xf[bhc_mask&s0_mask]['Gene Symbol']
# xf.loc[bhc_mask&s0_mask&bhc_d_mask, f'{cat_code} Sig FDR'] = 'decreased expression'

kmask = df['Gene Symbol'].str.upper().isin(known)
pmask = df['Gene Symbol'].str.upper().isin(prediction)

df['e3 substrate'] = ''
df['e3label'] = ''
df['k-e3 substrate'] = ''
df['k-e3label'] = ''


df.loc[kmask, 'e3 substrate'] = 'known substrate'
df.loc[pmask, 'e3 substrate'] = 'predicted substrate'
df.loc[kmask, 'e3label'] = df[kmask]['Gene Symbol']
df.loc[pmask, 'e3label'] = df[pmask]['Gene Symbol']
df.loc[kmask, 'k-e3 substrate'] = 'known substrate'
df.loc[kmask, 'k-e3label'] = df[kmask]['Gene Symbol']

df.to_csv(filepath+file.split('.')[0]+'_itxomeK.csv', index=False)




























