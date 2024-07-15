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

membrane_atlas = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\13. membrane proteomics\\'+
                  'NatCancer2021_surfaceProteome_tableS2.csv')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['Gene'].values)


path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\SGCellSurface\\'
match = 'Accession'


print('Reading df1: ', end='')
dfC = pd.read_excel(path+'sudhl8 DIA.xlsx')
print(f'{dfC.shape[0]} rows. {dfC.shape[1]} columns.')
print('Reading df2: ', end='')
dfS = pd.read_excel(path+'sudhl8 procan.xlsx')
print(f'{dfS.shape[0]} rows. {dfS.shape[1]} columns.')
print()

# print('Isolating Surface Proteins:')
# dfS = dfS[dfS['Gene Symbol'].isin(matlas)]
# dfC = dfC[dfC['Gene Symbol'].isin(matlas)]
# print(f'df1: {dfC.shape[0]} rows. {dfC.shape[1]} columns.')
# print(f'df2: {dfS.shape[0]} rows. {dfS.shape[1]} columns.')


print('Log2Transforming df1')
dfCdata = dfC[dfC.columns[2:]].replace(0, np.NaN)
dfCdata = np.log2(dfCdata/(dfCdata.median()))
# print('Log2Transforming df2')
# dfSdata = dfS[dfS.columns[3:]].replace(0, np.NaN)
# dfSdata = np.log2(dfSdata/(dfSdata.median()))

# for col in dfSdata.columns:
#     dfS[col] = dfSdata[col]

for col in dfCdata.columns:
    dfC[col] = dfCdata[col]


print()
mintersect = [x for x in set(dfS[match].unique()).intersection(set(dfC[match].unique()))]
dfS = dfS[dfS[match].isin(mintersect)]
dfC = dfC[dfC[match].isin(mintersect)]
print(f'Selecting intersecting proteins: {len(mintersect)} unique genes')
print(f'df1: {dfC.shape[0]} rows. {dfC.shape[1]} columns.')
print(f'df2: {dfS.shape[0]} rows. {dfS.shape[1]} columns.')
print()


dfS = dfS.sort_values(by=[match])
dfC = dfC.sort_values(by=[match])

# dfS.to_excel(path+'SGCellSurface\\exports\\dfS--2.xlsx',index=False)
# dfC.to_excel(path+'SGCellSurface\\exports\\dfC--2.xlsx',index=False)
# dfS = pd.read_excel(path+'SGCellSurface\\exports\\dfS--2.xlsx')
# dfC = pd.read_excel(path+'SGCellSurface\\exports\\dfC--2.xlsx')

cl_intersect = [x for x in set(dfS.columns[3:]).intersection(set(dfC.columns[2:]))]
print(f"Selecting intersecting cell lines: {len(cl_intersect)} unique lines")


xS = []
yC = []
clmen = []
prmen = []


print('Calculating matching Protein:Cell.Line pairs')
for cl in cl_intersect:
    sliceS = pd.Series([x for x in dfS[cl]], index=dfS[match])
    sliceC = pd.Series([x for x in dfC[cl]], index=dfC[match])

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

maSCdf = pd.DataFrame({'Cell Line':clmen, match:prmen, 'df1':yC, 'df2':xS})
maSCdf = maSCdf.dropna(subset=['df1', 'df2', match])

# maSCdf.to_excel(path+'maSCdf--.xlsx',index=False)

length = maSCdf.shape[0]
x = maSCdf['df1'].values.reshape(length, 1)
y = maSCdf['df2'].values.reshape(length, 1)
regr = linear_model.LinearRegression()
regr.fit(x, y)

score = r2_score(y, regr.predict(x))

plt.title(f"lin r2: {round(score, 4)}")
plt.ylabel('SUDHL8 procan')
plt.xlabel('SUDHL8 LAJ DIA')
plt.scatter(x, y, color='blue', s=1)
plt.plot(x, regr.predict(x), color='black', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()































