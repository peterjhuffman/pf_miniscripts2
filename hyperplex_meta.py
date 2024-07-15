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

inf = 1200
code = 'WT'
start = time()
filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\19. hyperplex meta\\'

#df_tmt = pd.read_excel(filepath+'PIK3CA_WT_TMT.xlsx')
#df_silac = pd.read_excel(filepath+'PIK3CA_WT_SILAC.xlsx')
#
#print(df_tmt.shape)
#print(df_silac.shape)

hdf = pd.read_csv(filepath+'WT_protein_intensity_normalized.csv')

# SINGLE PROTEIN - VISUAL


ACC = 'O95785'

xhdf = hdf[hdf['ProteinID']==ACC]


h_xhdf = xhdf[xhdf['channel']=='heavy']
l_xhdf = xhdf[xhdf['channel']=='light']


xh_inf = h_xhdf.iloc[0, :]
xl_inf = l_xhdf.iloc[0, :]
xh_inf['value'],xl_inf['value']  = 1,0
xh_inf['time'],xl_inf['time'] = inf,inf
for _ in range(3):
    h_xhdf = h_xhdf.append(xh_inf)
    l_xhdf = l_xhdf.append(xl_inf)

fig, (ax1, ax2) = plt.subplots(1,2,sharey=True)

fig.suptitle(f'Half Life of {xhdf["Gene"].values[0]} : {ACC}')


ax1.scatter(l_xhdf['time'], l_xhdf['value'], label='light')
ax1.scatter(h_xhdf['time'], h_xhdf['value'], label='heavy')
ax1.legend()


ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax1.set_xticks([x for x in h_xhdf['time'].unique()])
ax2.set_xticks([x for x in h_xhdf['time'].unique()])

ax1.set_ylim(-0.05,1.05)
ax2.set_ylim(-0.05,1.05)
ax2.set_xlim(-2,28)
ax1.set_xlim(-2,28)

#ax1.title(f'Half Life of {xhdf["Gene"].values[0]} : {ACC}')
#ax1.set_ylabel('heavy/light ratio')
#ax1.set_xlabel('time (h)')
#plt.show()

h_p0 = [-0.5, 0, 0.8]
l_p0 = [0.5, 0, 0.2]
x = np.linspace(0,100,200)

h_popt, h_pcov = curve_fit(func, h_xhdf['time'], h_xhdf['value'], p0=h_p0)
l_popt, l_pcov = curve_fit(func, l_xhdf['time'], l_xhdf['value'], p0=l_p0)
ax2.plot(x, func(x, *l_popt))
ax2.plot(x, func(x, *h_popt))
ax2.scatter(l_xhdf['time'], l_xhdf['value'], label='light')
ax2.scatter(h_xhdf['time'], h_xhdf['value'], label='heavy')
#ax2.ylabel('heavy/light ratio')
#ax2.ylabel('time (h)')

hl = 0-np.log((0.5 - h_popt[2])/h_popt[0])/h_popt[1]
ax2.axvline(x = hl, color = 'g')
ax1.set(xlabel='time (h)', ylabel='heavy/light ratio')
ax2.set(xlabel=f'Half-life: {round(hl,1)} hr')
plt.savefig(filepath+f'pulls\\wiz.png',dpi=400)
#print(f'Saved to {filepath}')
plt.show()


# BATCH - COMPUTE RAW
#
#
#ACCs = hdf['ProteinID'].unique()#[:500]
#GENEs = []
#hls = []
#hlacc = []
#
#num = 0
#denom = len(ACCs)
#for ACC in ACCs:
#    num += 1
#    print(f"{num}/{denom} :",end='')
#    xhdf = hdf[hdf['ProteinID']==ACC]
#    GENEs.append(xhdf['Gene'].values[0])
#
#    h_xhdf = xhdf[xhdf['channel']=='heavy']
#    l_xhdf = xhdf[xhdf['channel']=='light']
#
#
#    try:
#        h_popt, h_pcov = curve_fit(func, h_xhdf['time'], h_xhdf['value'])
#        l_popt, l_pcov = curve_fit(func, l_xhdf['time'], l_xhdf['value'])
#        hl = 0-np.log((0.5 - h_popt[2])/h_popt[0])/h_popt[1]
#        print('o')
#        if np.diag(h_pcov).sum()<0:
#            raise RuntimeError('terrible fit')
#        hl_meta = np.sqrt(np.diag(h_pcov)).sum()
#    
#    except RuntimeError:
#        hl = 'no convergence'
#        hl_meta = 'no convergence'
#        print('x')
#    
#    hls.append(hl)
#    hlacc.append(hl_meta)
#
#
#finF = pd.DataFrame({'Accession':ACCs, 'Gene':GENEs, 'Half-life':hls, 'Half-life estcov':hlacc})
#finF.to_excel(filepath+f'{code}_hlSILAC_cut_raw.xlsx',index=False)
#del finF
#
#
#print(f"time: {time()-start}")
#



# BATCH - COMPUTE ASSUME d_inf=1
#
#
#ACCs = hdf['ProteinID'].unique()#[:500]
#GENEs = []
#hls = []
#hlacc = []
#
#num = 0
#denom = len(ACCs)
#for ACC in ACCs:
#    num += 1
#    print(f"{num}/{denom} :",end='')
#    xhdf = hdf[hdf['ProteinID']==ACC]
#    GENEs.append(xhdf['Gene'].values[0])
#
#    h_xhdf = xhdf[xhdf['channel']=='heavy']
#    l_xhdf = xhdf[xhdf['channel']=='light']
#
#    xh_inf = h_xhdf.iloc[0, :]
#    xl_inf = l_xhdf.iloc[0, :]
#    xh_inf['value'],xl_inf['value']  = 1,0
#    xh_inf['time'],xl_inf['time'] = inf,inf
#    for _ in range(3):
#        h_xhdf = h_xhdf.append(xh_inf)
#        l_xhdf = l_xhdf.append(xl_inf)
#
#    try:
#        h_popt, h_pcov = curve_fit(func, h_xhdf['time'], h_xhdf['value'])
#        l_popt, l_pcov = curve_fit(func, l_xhdf['time'], l_xhdf['value'])
#        hl = 0-np.log((0.5 - h_popt[2])/h_popt[0])/h_popt[1]
#        print('o')
#        if np.diag(h_pcov).sum()<0:
#            raise RuntimeError('terrible fit')
#        hl_meta = np.sqrt(np.diag(h_pcov)).sum()
#
#    except RuntimeError:
#        hl = 'no convergence'
#        hl_meta = 'no convergence'
#        print('x')
#
#    hls.append(hl)
#    hlacc.append(hl_meta)
#
#
#finF = pd.DataFrame({'Accession':ACCs, 'Gene':GENEs, 'Half-life':hls, 'Half-life estcov':hlacc})
#finF.to_excel(filepath+f'{code}_hlSILAC_asm.xlsx',index=False)
#del finF
#
#
#print(f"time: {time()-start}")
#
#
#
#
### BATCH - COMPUTE VARIABLE ASSUMPTION
#
#
#ACCs = hdf['ProteinID'].unique()#[:500]
#GENEs = []
#hls = []
#hlacc = []
#
#num = 0
#denom = len(ACCs)
#for ACC in ACCs:
#    num += 1
#    print(f"{num}/{denom} :",end='')
#    xhdf = hdf[hdf['ProteinID']==ACC]
#    GENEs.append(xhdf['Gene'].values[0])
#
#    h_xhdf = xhdf[xhdf['channel']=='heavy']
#    l_xhdf = xhdf[xhdf['channel']=='light']
#
#
#    try:
#        h_popt, h_pcov = curve_fit(func, h_xhdf['time'], h_xhdf['value'])
#        l_popt, l_pcov = curve_fit(func, l_xhdf['time'], l_xhdf['value'])
#        hl = 0-np.log((0.5 - h_popt[2])/h_popt[0])/h_popt[1]
#        print('o')
#        if np.diag(h_pcov).sum()<0:
#            raise RuntimeError('terrible fit')
#        hl_meta = np.sqrt(np.diag(h_pcov)).sum()
#    
#    except RuntimeError:
#        try:
#            xh_inf = h_xhdf.iloc[0, :]
#            xl_inf = l_xhdf.iloc[0, :]
#            xh_inf['value'],xl_inf['value']  = 1,0
#            xh_inf['time'],xl_inf['time'] = inf,inf
#            for _ in range(3):
#                h_xhdf = h_xhdf.append(xh_inf)
#                l_xhdf = l_xhdf.append(xl_inf)
#                h_popt, h_pcov = curve_fit(func, h_xhdf['time'], h_xhdf['value'])
#                l_popt, l_pcov = curve_fit(func, l_xhdf['time'], l_xhdf['value'])
#                hl = 0-np.log((0.5 - h_popt[2])/h_popt[0])/h_popt[1]
##                print('o')
#                if np.diag(h_pcov).sum()<0:
#                    raise RuntimeError('terrible fit')
#                hl_meta = np.sqrt(np.diag(h_pcov)).sum()
#        except RuntimeError:
#            hl = 'no convergence'
#            hl_meta = 'no convergence'
#            print('x')
#
#    hls.append(hl)
#    hlacc.append(hl_meta)
#
#
#finF = pd.DataFrame({'Accession':ACCs, 'Gene':GENEs, 'Half-life':hls, 'Half-life estcov':hlacc})
#finF.to_excel(filepath+f'{code}_hlSILAC_var.xlsx',index=False)
#del finF
#
#
#print(f"time: {time()-start}")
#
#










