# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 09:55:27 2024

@author: HUFFMP
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df22 = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\kat2global\\exports\\PFEC_PH_20231026_K2Lglobal_fx_prot.xlsx')
membrane_atlas = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\membrane proteomics\\'+
                  'NatCancer2021_surfaceProteome_tableS2.csv')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['Gene'].values)

#dfSurface = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\SGCellSurface\\SGCellSurface_surface.xlsx')
#dfTurno = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\SGCellSurface\\SGCellSurface_turnover.xlsx')

df22abcol = [x for x in df22.columns if str(x).count('Abundance:')>0]
df22['Abundance Total'] = df22[df22abcol].sum(axis=1)
df22['Log Abundance Total'] = np.log2(df22['Abundance Total'])



df22surf = df22[df22['Gene Symbol'].isin(matlas)]

print(df22surf.shape)



df22surfsort = df22surf.sort_values(by=['Abundance Total'], ascending=False)
df22surfsort['ab Rank'] = [x+1 for x in range(df22surfsort.shape[0])]

plt.scatter(df22surfsort['ab Rank'], df22surfsort['Log Abundance Total'])
plt.suptitle('Log2 Rank Plot of Membrane Protein Abundance in 22rv1 Global')
plt.title(f'22rv1:{df22surf.shape[0]} srf')
plt.ylabel('Log2 Total Protein Abundance')
plt.xlabel('Rank Protein Abundance')
plt.show()


cell_target = 'LNCAP'

sea_target = dfSurface[dfSurface[cell_target]!=0]
prot_target = [x for x in sea_target['GeneSymbol']]
target_seasrf = sea_target[sea_target['GeneSymbol'].isin(matlas)].shape[0]


df22sgN = df22surfsort[~df22surfsort['Gene Symbol'].isin(prot_target)]
df22sgY = df22surfsort[df22surfsort['Gene Symbol'].isin(prot_target)]

plt.figure(figsize=(25,5))
plt.scatter(df22sgN['ab Rank'], df22sgN['Log Abundance Total'], s=1, marker='o', label='SG miss')
plt.scatter(df22sgY['ab Rank'], df22sgY['Log Abundance Total'], s=1, marker='o', label='SG hit')
plt.suptitle('Log2 Rank Plot of Membrane Protein Abundance in LNCAP Global')
plt.title(f'LNCAP:{df22surf.shape[0]} srf, SEA-{cell_target}:{target_seasrf} srf')
plt.ylabel('Log2 Total Protein Abundance')
plt.xlabel('Rank Protein Abundance')
plt.legend()
plt.savefig(f'srf_rank_{cell_target}_ln.png', dpi=400)
plt.show()



#df22sgY.to_excel('df22sgY.xlsx')











































