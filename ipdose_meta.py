# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:33:45 2023

@author: HUFFMP
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from time import time
pd.options.mode.chained_assignment = None  # default='warn'

def func(x, a, b, c):
    return a * np.exp(-b * x) + c


def rep2_interp(df):
    itercol = df.columns
    df['interpolated col'] = '-'
    
    switch = True
    for val, col in enumerate(itercol):
        if switch:
            copycol = df.columns[val+1]
        else:
            copycol = df.columns[val-1]
        df.loc[df[col]==0.0, 'interpolated col'] = df.loc[df[col]==0.0, 'interpolated col'].str.replace('-', f'{val+1}, -')
        df.loc[df[col]==0.0, col] = df.loc[df[col]==0.0, copycol]
        switch = switch==False
    return df

def hill(x, a, b, c):
    y = a + (c*(x/(b+x)))
    return y

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(k*(x+x0))) + b
    return (y)



filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\jzd\\'


df = pd.read_excel(filepath+'8951_5211_MCF7_Depletion_Factor_PDoutput.xlsx')

trash = ['Abundance: F1: 126, Sample','Abundance: F1: 127N, Sample', 'Abundance: F1: 128C, Sample',
         'Abundance: F1: 129N, Sample','Abundance: F1: 130C, Sample', 'Abundance: F1: 131N, Sample',
         'Abundance: F1: 132C, Sample', 'Abundance: F1: 133N, Sample','Abundance: F1: 134C, Sample', 'Abundance: F1: 135N, Sample']
for x in trash:
    del df[x]


abcols = [f for f in df.columns if f.count('Abundance:')>0]

# dropping unquant rows
hdf = df.dropna(subset = abcols)
print(hdf.shape)

hdf['top ab'] = hdf[abcols].max(axis=1)
AB_CUTOFF = 1000

hdf = hdf[hdf['top ab']>=AB_CUTOFF]
del hdf['top ab']
print(hdf.shape)

else_cols = list(set(hdf.columns).difference(set(abcols)))
df_aux = hdf[else_cols]

hdf_med = hdf[abcols].median()
hdf_main = hdf[abcols]/hdf_med
print(hdf_main.columns)

#hdf_main.to_csv(filepath+'input.csv', index=False)
#hdf_main = rep2_interp(hdf_main)
#hdf_main.to_csv(filepath+'output.csv', index=False)
#df_aux['interpolated col'] = hdf_main['interpolated col']
#del hdf_main['interpolated col']

limit_of_detection = hdf_main.replace(0, 100).min()[abcols].min()
hdf_main = hdf_main.replace(0, limit_of_detection)


hdf_tf = np.log2(hdf_main)
hdf_tf = hdf_tf-hdf_tf.min().min()
#hdf_tf = hdf_main

for col in else_cols:
    hdf_tf[col] = df_aux[col]

doses = ['10uM', '1uM', '100nM', 'DMSO']
hdf_tf[doses[0]] = (hdf_tf['Abundance: F1: 129C, Sample']+hdf_tf['Abundance: F1: 130N, Sample'])/2
hdf_tf[doses[1]] = (hdf_tf['Abundance: F1: 131C, Sample']+hdf_tf['Abundance: F1: 132N, Sample'])/2
hdf_tf[doses[2]] = (hdf_tf['Abundance: F1: 133C, Sample']+hdf_tf['Abundance: F1: 134N, Sample'])/2
hdf_tf[doses[3]] = (hdf_tf['Abundance: F1: 127C, Sample']+hdf_tf['Abundance: F1: 128N, Sample'])/2

#semi = True
GENEs = ['Q16659']
#GENEs = [x for x in hdf_tf['Accession']]

notabs = [3,2,1,0]
tabs = [10, 1, 0.1, 0]
logtab = [np.log10(x+0.00001) for x in tabs]

oc50s = []
pcovs = []

for target in GENEs:
    print('-', end='')
    gvals = hdf_tf.loc[hdf['Accession']==target, doses].values[0]
#    gvals = [2**x for x in gvals]
#    print(gvals)
#    plt.scatter(doses, gvals)
#    plt.title(target)
#    plt.ylabel('abundance logtf')
#    plt.xlabel('lane (not blocked)')
#    plt.show()
#    plt.clf()
    plt.scatter(tabs, gvals)
    plt.ylabel('abundance logtf')
    plt.xlabel('lane (not blocked)')
#    plt.show()
#    plt.clf()
#    plt.scatter([np.log2(x+0.000001) for x in [10, 1, 0.1, 0]], gvals)
#    plt.title(target)
#    plt.ylabel('abundance logtf')
#    plt.xlabel('lane (not blocked)')
#    plt.show()
#    plt.clf()

    try:
        h_popt, h_pcov = curve_fit(func, tabs, gvals)


    #, p0=h_p0
        x = np.linspace(0,15,200)
    #    hl_y = gvals[0] + (gvals[3]-gvals[0])/2
        hl_y = gvals[3]-1
#        hl =1
        hl = np.log((hl_y-h_popt[2])/h_popt[0]) / -h_popt[1]
        corr = round(np.diag(h_pcov).sum(),3)
        
        plt.plot(x, func(x, *h_popt))
        plt.title(f"{target}     oc50: {int(round(hl*1000, 0))}nM    pcov:{corr}")
        plt.show()
        plt.clf()
        if (hl>0) and (corr<30):
            oc50s.append(hl)
            pcovs.append(np.diag(h_pcov).sum())
        else:
            oc50s.append('no fit')
            pcovs.append('no fit')
    except:
        oc50s.append('no fit')
        pcovs.append('no fit')

#hdf_tf['OC50'] = oc50s
#hdf_tf['pcovariance'] = pcovs
#
#hdf_tf.to_excel(filepath+'hls_modeled.xlsx', index=False)



