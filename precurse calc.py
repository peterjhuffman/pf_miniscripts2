# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:12:18 2024

@author: HUFFMP
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit





def list_files_in_directory(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def list_folders(directory_path):
    folders = []
    for item in os.listdir(directory_path):
        item_path = os.path.join(directory_path, item)
        if os.path.isdir(item_path):
            folders.append(item)
    return folders

def folderslicer(fraws):
    fraws = [f for f in fraws if len(f.split('_'))>4]
    fraws = [f for f in fraws if str(f).lower().count('wash')==0]
    
    numcaps = list(set(['_'.join(f.split('_')[:min(4, len(f.split('_')))]) for f in fraws]))

    foldproj = []
    numcounter = [len([raw for raw in fraws if raw.count(cap)>0]) for cap in numcaps]

    for cap,count in zip(numcaps,numcounter):
        if count>3:
            stocklimit=4
            bencher = fraws[[f.count(cap) for f in fraws].index(1)]
            while len([raw for raw in fraws if raw.count('_'.join(bencher.split('_')[:min(stocklimit, len(bencher.split('_')))]))>0])==count:
                stocklimit+=1
            foldproj.append('_'.join(bencher.split('_')[:min(stocklimit-1, len(bencher.split('_')))]))
    return foldproj


def file_in_dir(file, directory):
    for root,dirs,files in os.walk(directory):
        for name in files:
            if name == file:
                return True
    return False

def gaussianfx(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


path = '//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/PFAS_AHH_240711_KIF18global_Fr1ms.xlsx'


scores = []
gaussdf = pd.DataFrame()
bincount = 80
rankr1 = 2.6
featurenum = 7
ilen = 6
savind=6

print('===>',end='')
# fold = 'PFAS_YYL_20240329_1INSM1_f7'
df = pd.read_excel(path)
print(path)
fold = '.'.join(path.split('.')[:-1])
plt.hist(df['Precursor Intensity'], bins=100, range=(0,600000))
plt.ylim([0, 2000])
plt.title(f"{fold} Precursor Intensities")
plt.ylabel('# Precursors')
plt.xlabel('Intensity')
# plt.savefig(f'{path}/precurse/{fold}_pci.png', dpi=300)
# plt.show()
plt.clf()

vals = plt.hist(df['Precursor Intensity'], bins=bincount, range=(0,600000))
plt.clf()
run = vals[1][1]-vals[1][0]

stock = vals[0]
s1 = np.sum(stock)


miniscore = [i*i*j for i,j in zip([x for x in range(len(vals[0]))], [x for x in vals[0]])]

print(f"PreCurse: {np.sum(miniscore/s1)}")
