# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:58:33 2024

@author: HUFFMP
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from sklearn import datasets, linear_model
from sklearn.metrics import r2_score
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

membrane_atlas = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\membrane proteomics\\'+
                  'NatCancer2021_surfaceProteome_tableS2.csv')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['Gene'].values)


path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\'

print('Reading CCLE: ', end='')
dfC = pd.read_excel(path+'ccle pandar.xlsx')
# dfC = pd.read_excel(path+'SGCellSurface\\PFEC_PH_20231026_K22global_controlC.xlsx')
print(f'{dfC.shape[0]} rows. {dfC.shape[1]} columns.')
print('Reading SGSurface: ', end='')
dfS = pd.read_excel(path+'SGCellSurface\\procan hcc1806.xlsx')
# dfS = pd.read_excel(path+'SGCellSurface\\SH_22rv1.xlsx')
print(f'{dfS.shape[0]} rows. {dfS.shape[1]} columns.')
print()

print('Isolating Surface Proteins:')
dfS = dfS[dfS['GeneSymbol'].isin(matlas)]
dfC = dfC[dfC['Gene Symbol'].isin(matlas)]
print(f'CCLE: {dfC.shape[0]} rows. {dfC.shape[1]} columns.')
print(f'SGSurface: {dfS.shape[0]} rows. {dfS.shape[1]} columns.')

print()
print('Log2Transforming SGSurface Data')
dfSdata = dfS[dfS.columns[3:]].replace(0, np.NaN)
dfSdata = np.log2(dfSdata/(dfSdata.median()))

dfCtrash = [x for x in dfC.columns[6:48]]
for trash in dfC.columns[6:48]:
    del dfC[trash]

dfC.columns = [str(x).split('_')[0] for x in dfC.columns]

for col in dfSdata.columns:
    dfS[col+' _px'] = dfSdata[col]


print()
mintersect = [x for x in set(dfS['GeneSymbol'].unique()).intersection(set(dfC['Gene Symbol'].unique()))]
dfS = dfS[dfS['GeneSymbol'].isin(mintersect)]
dfC = dfC[dfC['Gene Symbol'].isin(mintersect)]
print(f'Selecting intersecting proteins: {len(mintersect)} unique genes')
print(f'CCLE: {dfC.shape[0]} rows. {dfC.shape[1]} columns.')
print(f'SGSurface: {dfS.shape[0]} rows. {dfS.shape[1]} columns.')
print()


dfS = dfS.sort_values(by=['GeneSymbol'])
dfC = dfC.sort_values(by=['Gene Symbol'])

# dfS.to_excel(path+'SGCellSurface\\exports\\dfS--2.xlsx',index=False)
# dfC.to_excel(path+'SGCellSurface\\exports\\dfC--2.xlsx',index=False)
# dfS = pd.read_excel(path+'SGCellSurface\\exports\\dfS--2.xlsx')
# dfC = pd.read_excel(path+'SGCellSurface\\exports\\dfC--2.xlsx')

cl_intersect = [x for x in set(dfS.columns[3:]).intersection(set(dfC.columns[6:]))]
print(f"Selecting intersecting cell lines: {len(cl_intersect)} unique lines")


xS = []
yC = []
clmen = []
prmen = []


print('Calculating matching Protein:Cell.Line pairs')
for cl in cl_intersect:
    sliceS = pd.Series([x for x in dfS[cl+' _px']], index=dfS['GeneSymbol'])
    sliceC = pd.Series([x for x in dfC[cl]], index=dfC['Gene Symbol'])

    # print(cl)

    for pr, y in zip(sliceC.index, sliceC):
        yC.append(y)
        if type(sliceS[pr]) == np.float64:
            xS.append(sliceS[pr])
        else:
            xS.append(sliceS[pr].max())
        # print(type(sliceS[pr]))
        clmen.append(cl)
        prmen.append(pr)

print('calculation complete.')

maSCdf = pd.DataFrame({'Cell Line':clmen, 'Gene Symbol':prmen, 'CCLE score':yC, 'SGSurface score':xS})
maSCdf = maSCdf.dropna(subset=['CCLE score', 'SGSurface score', 'Gene Symbol'])

maSCdf.to_excel(path+'SGCellSurface\\exports\\maSCdf--.xlsx',index=False)

length = maSCdf.shape[0]
x = maSCdf['CCLE score'].values.reshape(length, 1)
y = maSCdf['SGSurface score'].values.reshape(length, 1)
regr = linear_model.LinearRegression()
regr.fit(x, y)

score = r2_score(y, regr.predict(x))

plt.title(f"r2: {round(score, 4)}")
plt.ylabel('SH 9xplex A375')
plt.xlabel('CCLE match')
plt.scatter(x, y, color='blue', s=1)
plt.plot(x, regr.predict(x), color='black', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()

































