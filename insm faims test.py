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
from time import sleep as slp



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

path = '//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/hela faims'
membrane_atlas = ('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/results/13. membrane proteomics/NatCancer2021_surfaceProteome_tableS2_illustrated.csv')
dfmatlas = pd.read_csv(membrane_atlas)
matlas = list(dfmatlas['accession'].values)
# slp(60*12)
prot = []
pep = []
psm = []
memprot = []
for fold in list_folders(path):
    print(fold)
    df = pd.read_excel(path+f'//{fold}//{fold}_prot.xlsx')
    prot.append(df[df['Master']=='IsMasterProtein'].shape[0])
    # df = df[df['Master']=='IsMasterProtein']
    # memprot.append(df[df['Accession'].isin(matlas)].shape[0])
    df = pd.read_excel(path+f'//{fold}//{fold}_pep.xlsx')
    pep.append(df.shape[0])
    df = pd.read_excel(path+f'//{fold}//{fold}_psm.xlsx')
    psm.append(df.shape[0])
    del df

gf = pd.DataFrame({'Run': list_folders(path), 'Proteins':prot, 'Peptides':pep, 'PSMs':psm})
gf.to_excel(f'{path}/global_scorecard.xlsx')


# df = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/PFAS_YYL_20240329_1INSM1_901_f1/PFAS_YYL_20240329_1INSM1_901_f1_msms.xlsx')



# for fold in list_folders(path):
#     print('===>',end='')
#     df = pd.read_excel(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/{fold}/{fold}_msms.xlsx')
#     print(fold)
#     plt.hist(df['Precursor Intensity'], bins=100, range=(0,600000))
#     plt.ylim([0, 2000])
#     plt.title(f"{fold} Precursor Intensities")
#     plt.ylabel('# Precursors')
#     plt.xlabel('Intensity')
#     # plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/precurse/{fold}_pci.png', dpi=300)
#     plt.show()
#     plt.clf()

#     plt.hist(df['Isolation Interference in Percent'], bins=100, range=(0,100))
#     plt.ylim([0, 5000])
#     plt.title(f"{fold} Isolation Interference")
#     plt.ylabel('# MSMS')
#     plt.xlabel('Intensity')
#     plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/postcurse/{fold}_ii.png', dpi=300)
#     # plt.show()
#     plt.clf()





# fold = 'PFEC_PH_20240624_ISDglobal_fx'
# df = pd.read_excel(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/{fold}/{fold}_msms.xlsx')
# print(fold)
# plt.hist(df['Precursor Intensity'], bins=100, range=(0,600000))
# plt.ylim([0, 2000])
# plt.title(f"{fold} Precursor Intensities")
# plt.ylabel('# Precursors')
# plt.xlabel('Intensity')
# plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/{fold}/{fold}_pci.png', dpi=300)
# # plt.show()
# plt.clf()

# scores = []
# gaussdf = pd.DataFrame()
# bincount = 80
# rankr1 = 2.6
# featurenum = 7
# ilen = 6
# savind=6
# foldcache = list_folders(path)[1:-2]

# for fold in foldcache:
#     print('===>',end='')
# # fold = 'PFAS_YYL_20240329_1INSM1_f7'
#     df = pd.read_excel(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/{fold}/{fold}_msms.xlsx')
#     print(fold)
#     # plt.hist(df['Precursor Intensity'], bins=100, range=(0,600000))
#     # plt.ylim([0, 2000])
#     # plt.title(f"{fold} Precursor Intensities")
#     # plt.ylabel('# Precursors')
#     # plt.xlabel('Intensity')
#     # # plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/precurse/{fold}_pci.png', dpi=300)
#     # # plt.show()
#     # plt.clf()
    
#     vals = plt.hist(df['Precursor Intensity'], bins=bincount, range=(0,600000))
#     plt.clf()
#     run = vals[1][1]-vals[1][0]
    
    
#     x = np.linspace(-1, 1, 50)
#     gaussian = norm.pdf(x, 0, 1)
    
#     stock = vals[0]
#     s1 = np.sum(stock)
    
#     rankr = (rankr1, 1)
#     X1s = []
#     Y1s = []
#     sig1S = []
    
#     for j in range(featurenum):
#         corr = signal.correlate(stock, gaussian, mode='same')
        
#         X1 = np.argmax((corr*rankr[1])+([s*rankr[0] for s in stock]))
    
#         if X1 in X1s:
#             X1 = max(np.argmax(stock), X1-5)
        
#         if stock[X1] == 0:
#             X1 = max(np.argmax(stock), X1-5)
        
#         Y1 = stock[X1]
        
        
#         sigmas = []
        
#         i = 1
        
#         while (i<ilen) and (X1+i<len(vals[0])):
#             X2 = X1 + i
#             Y2 = stock[X1+i]
#             sigmas.append(np.sqrt(-0.5 * (X2 - X1)**2 / np.log(Y2 / Y1)))
#             i+=1
        
#         sigma_est=np.mean([x for x in pd.Series(sigmas).dropna().values])

#         b = 0
#         ilenbox = [2,4,5,7,9,11,13]
#         while (np.isnan(sigma_est) or sigma_est>12) and b<7:
#             ilen = ilenbox[b]
#             # print(ilen)
#             while (i<(ilen)) and (X1+i<len(vals[0])):
#                 X2 = X1 + i
#                 Y2 = stock[X1+i]
#                 sigmas.append(np.sqrt(-0.5 * (X2 - X1)**2 / np.log(Y2 / Y1)))
#                 i+=1
            
#             sigma_est=np.mean([x for x in pd.Series(sigmas).dropna().values])
#             b +=1

#         bi = 0.1
#         while np.isnan(sigma_est) or sigma_est>12:
#             print(bi)
#             ranked = 0+bi
#             X1 = np.argmax((corr*rankr[1])+([s*ranked for s in stock]))
        
#             if X1 in X1s:
#                 X1 = max(np.argmax(stock), X1-5)
            
#             if stock[X1] == 0:
#                 X1 = max(np.argmax(stock), X1-5)
            
#             Y1 = stock[X1]
            
            
#             sigmas = []
            
#             i = 1
            
#             while (i<(ilen+(X1/14))) and (X1+i<len(vals[0])):
#                 X2 = X1 + i
#                 Y2 = stock[X1+i]
#                 sigmas.append(np.sqrt(-0.5 * (X2 - X1)**2 / np.log(Y2 / Y1)))
#                 i+=1
            
#             sigma_est=np.mean([x for x in pd.Series(sigmas).dropna().values])

#             b = 0
#             ilenbox = [2,4,5,7,9,11,13]
#             while (np.isnan(sigma_est) or sigma_est>12) and b<7:
#                 ilen = ilenbox[b]
#                 # print(ilen)
#                 while (i<(ilen+(X1/12))) and (X1+i<len(vals[0])):
#                     X2 = X1 + i
#                     Y2 = stock[X1+i]
#                     sigmas.append(np.sqrt(-0.5 * (X2 - X1)**2 / np.log(Y2 / Y1)))
#                     i+=1
                
#                 sigma_est=np.mean([x for x in pd.Series(sigmas).dropna().values])
#                 b +=1
#             bi += 0.1

#         X1s.append(X1)
#         Y1s.append(Y1)
#         sig1S.append(sigma_est)
        
#         x = np.linspace(0, bincount, bincount)
#         y = gaussianfx(x, Y1, X1, sigma_est)
#         plt.plot(x, y)
#         plt.plot(stock)
#         plt.plot((corr/corr.max())*np.max(stock))
#         plt.title(f"#{j} {fold}")
#         plt.show()
#         plt.clf()

#         stock = [max(0, x-i) for x,i in zip(stock, y)]


#     plt.plot(stock)
#     plt.show
#     s2 = (s1 - np.sum(stock))/s1

#     minigauss = pd.DataFrame()

#     minigauss[f"{fold} fX"] = X1s
#     minigauss[f"{fold} fY"] = Y1s
#     minigauss[f"{fold} fS: {round(s2, 2)}"] = sig1S

#     minigauss = minigauss.sort_values(by=f"{fold} fX")
#     minigauss = minigauss.reset_index(drop=True)

#     for c in minigauss.columns:
#         gaussdf[c] = minigauss[c]

#     miniscore = np.sum(minigauss[f"{fold} fX"]*minigauss[f"{fold} fX"]*minigauss[f"{fold} fY"]*minigauss[f"{fold} fS: {round(s2, 2)}"])
    

#     gaussdf.to_excel(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/faims_gauss_v{savind}.xlsx',index=False)
#     scores.append(miniscore)

# dfi = pd.DataFrame({'file':[x for x in foldcache], 'score':[np.sum(score) for score in scores]})
# dfi.to_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/PD export/insm1 faims/faims_gauss_scoring.xlsx', index=False)
