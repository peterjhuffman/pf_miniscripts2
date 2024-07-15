# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 09:36:00 2022

@author: HUFFMP
"""

import pandas as pd
import sys
from time import sleep as slp
from matplotlib import pyplot as plt

def stop():
    print('Exiting mixratio.exe. Thanks for dropping by.')
    slp(2)
    sys.exit()

datasource = input('DATA LOCATION (x to exit): ')
print('Loading file.')
if datasource.lower()=='x':
    stop()
fdtype = datasource[-3:]
if fdtype == 'csv':
    df = pd.read_csv(datasource)
elif fdtype == 'lsx':
    df = pd.read_excel(datasource)
else:
    print('filetype not recognized.')
    stop()
print('File read successfully.\n')

print(df.shape)
#df = df[df['Quan Info']!='No Quan Values']
#df = df[df['Quan Info']!='Excluded by Method']
print(df.shape)

sav = '.'.join(datasource.split('.')[:-1]).split('\\')[-1]
newpath = f'\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\tmt mixing'



THRESH = str(input('THRESHOLD (ex: enter 0.15 for the middle 70%)(x to exit): '))
if THRESH.lower()=='x':
    stop()
THRESH = float(THRESH)
if (THRESH > 1) or (THRESH<0):
    print('Incorrect threshold notation. Ask Peter.')
    stop()

ab_count = 1
for ab in df.columns:
    print(f"{ab_count}: {ab}")
    ab_count+=1
print()
print('Abundance columns detected.')
ab_select = input('Enter the start and end indices of the desired abundance columns for analysis, separated by a comma.\n'+
                  'Example: if columns 1-18 are desired, input should be: 1,18.\n'+
                  'Type ALL for all columns.\n\nIndices: ')

abundances = df.columns
if ab_select.upper() == 'ALL':
    abundances = df.columns
elif len(ab_select.split(','))!=2:
    print('Incorrect Formatting.')
    stop()
else:
    ab_indc = ab_select.split(',')
    abundances = abundances[int(ab_indc[0])-1:int(ab_indc[1])]

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

means = []
for col in abundances:
    df_slice = df[col]
    df_slice = df_slice[(df_slice>df_slice.quantile(THRESH))&(df_slice<df_slice.quantile(1-THRESH))]
    means.append(df_slice.mean())

plt.bar(abundances, means)
plt.ylabel('Mean Abundance')
plt.xlabel('TMT Lane')
plt.xticks([i for i in range(0, 18)], [str(i) for i in range(1, 19)])
plt.title('TMT Lanes: Mean Abundances')
plt.savefig(newpath+f'/{sav}_bar.png', dpi=1200)
print(f'barplot saved to folder ChemBio\MS_Analyses_PD\TMTtesting/{sav}')
plt.clf()

tmtMXR = pd.DataFrame()
tmtMXR['means'] = means

meanmin = max(means)
for val in means:
    if (val < meanmin) and (val*4 > max(means)):
        meanmin = val

tmtMXR['rel. ratio'] = tmtMXR['means'] / meanmin
tmtMXR['add. percent'] = 1 / tmtMXR['rel. ratio']

tmtMXR['add. percent'][tmtMXR['add. percent']>1] = 0
tmtMXR = tmtMXR.round({'add. percent':3, 'rel. ratio':4, 'means':4})

tmtMXR.columns = ['Abundance Means', 'Relative Abundance Ratio', 'New Mix Ratio']

tmtMXR.to_csv(newpath+f'/{sav}_values.csv', index=False)

plt.table(cellText=tmtMXR.values, colLabels=tmtMXR.columns, loc='center', cellLoc='center')
plt.axis('off')
plt.title('TMT mix ratio calculation', y=1.08)
plt.tight_layout()
plt.xticks([])
plt.yticks([])
print(newpath+f'/{sav}_summary{str(int(THRESH*100))}.png')
plt.savefig(newpath+f'/{sav}_summary{str(int(THRESH*100))}.png', dpi=1200)
print(f'TMT mix ratio table saved to folder ChemBio\MS_Analyses_PD\TMTtesting/{sav}')
plt.clf()