# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:58:16 2023

@author: HUFFMP
"""

import numpy as np
import pandas as pd


datasourcemax = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\membrane proteomics\\'+
                  'PFEC_PH_041423_XMP_APOBEC_prot.xlsx')

membrane_atlas = ('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/results/13. membrane proteomics/NatCancer2021_surfaceProteome_tableS2.csv')


# dfm = pd.read_excel(datasourcemax)


# df22r = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\kat2global\\exports\\PFAS_PH_20231103_K22global_fx_prot.xlsx')
#
bruh = 'Labelfree_DDA_hela_12fxns_chimerys_proteins_filtered'

dfSEAGEN = pd.read_excel(f'\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\{bruh}.xlsx')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['Gene'].values)





# print('22RV1 Global ISFX:')
# print(f"Proteins:                   {df22r.shape[0]}")
# print(f"mAtlas Surface Proteins:    {df22r[df22r['Gene Symbol'].isin(matlas)].shape[0]} ({round(100*df22r[df22r['Gene Symbol'].isin(matlas)].shape[0]/df22r.shape[0], 1)}%)")
# print(f"mAtlas Surface Peptides:    {df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides'].sum()} (avg: {round((df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides'].sum())/(df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides'].size), 2)}) (med: {df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides'].median()}) (1pep prot: {df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides'][df22r[df22r['Gene Symbol'].isin(matlas)]['# Peptides']==1].size})")
# print()
# print()

# print("Membrane Phase-Enriched Proteomics:")
# print(f"Proteins:                   {dfm.shape[0]}")
# print(f"mAtlas Surface Proteins:    {dfm[dfm['Gene Symbol'].isin(matlas)].shape[0]} ({round(100*dfm[dfm['Gene Symbol'].isin(matlas)].shape[0]/dfm.shape[0], 1)}%)")
# print(f"mAtlas Surface Peptides:    {dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides'].sum()} (avg: {round((dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides'].sum())/(dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides'].size), 2)}) (med: {dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides'].median()}) (1pep prot: {dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides'][dfm[dfm['Gene Symbol'].isin(matlas)]['# Peptides']==1].size})")

print("Labelfree DDA Hela 12fx:")
print(f"Proteins:                   {dfSEAGEN.shape[0]}")
print(f"mAtlas Surface Proteins:    {dfSEAGEN[dfSEAGEN['Gene Symbol'].isin(matlas)].shape[0]}")
#print(f"mAtlas Surface Peptides:    {dfSEAGEN[dfSEAGEN['Gene Symbol'].isin(matlas)]['# Peptides'].sum()}")



#
#setpot1 = set(list(dfpot1[dfpot1['Gene Symbol'].isin(matlas)]['Gene Symbol'].values))
#setpot2 = set(list(dfpot2[dfpot2['Gene Symbol'].isin(matlas)]['Gene Symbol'].values))
#setpot3 = set(list(dfpot3[dfpot3['Gene Symbol'].isin(matlas)]['Gene Symbol'].values))
#setpot4 = set(list(dfpot4[dfpot4['Gene Symbol'].isin(matlas)]['Gene Symbol'].values))
#setpot5 = set(list(dfpot5[dfpot5['Gene Symbol'].isin(matlas)]['Gene Symbol'].values))
#
#print('union: ', end=' ')
#print(len(setpot1.union(setpot2.union(setpot3.union(setpot4.union(setpot5))))))
