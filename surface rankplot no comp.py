# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 09:55:27 2024

@author: HUFFMP
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df22 = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\Labelfree_DDA_hela_24fxns_proteins_filtered 1_Hruler.xlsx')
membrane_atlas = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\membrane proteomics\\'+
                  'NatCancer2021_surfaceProteome_tableS2.csv')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['Gene'].values)

#dfSurface = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\SGCellSurface\\SGCellSurface_surface.xlsx')
#dfTurno = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\SGCellSurface\\SGCellSurface_turnover.xlsx')

df22abcol = [x for x in df22.columns if str(x).count('copyNum')>0]
df22['Abundance Total'] = df22[df22abcol].sum(axis=1)
df22['Log Abundance Total'] = np.log2(df22['Abundance Total'])





df22sort = df22.sort_values(by=['Abundance Total'], ascending=False)
df22sort['ab Rank'] = [x+1 for x in range(df22sort.shape[0])]

plt.scatter(df22sort['ab Rank'], df22sort['Log Abundance Total'])
plt.suptitle('Log2 Rank Plot of Membrane Protein Abundance in HeLa labelfree')
plt.title(f'HeLa:{df22.shape[0]}')
plt.ylabel('Log2 Total Protein Abundance')
plt.xlabel('Rank Protein Abundance')
plt.show()




df22sgN = df22sort
df22sgY = df22sort[df22sort['Gene Symbol'].isin(matlas)]
df22sgY = pd.concat([df22sgY.iloc[:20][:],df22sgY.iloc[-20:][:]])

df22sgY.to_excel('srfranker_leaderboard.xlsx',index=False)

plt.figure(figsize=(30,5))
plt.scatter(df22sgN['ab Rank'], df22sgN['Log Abundance Total'], s=1, marker='o', label='Internal')
plt.scatter(df22sgY['ab Rank'], df22sgY['Log Abundance Total'], s=1, marker='o', label='Surface')
plt.suptitle('Log2 Rank Plot of Membrane Protein Abundance in LNCAP Global')
plt.title(f'HeLa:{df22.shape[0]}, {df22sgY.shape[0]}srf')
plt.ylabel('Log2 Total Protein copynum')
plt.xlabel('Rank Protein Abundance')
plt.legend()
plt.savefig(f'srf_rank_hela_copyn.png', dpi=700)
plt.show()



#df22sgY.to_excel('df22sgY.xlsx')











































