# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:58:38 2024

@author: huffmp
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import re
import requests
from requests.adapters import HTTPAdapter, Retry
from scipy import stats
import math
from scipy.optimize import curve_fit


pd.options.mode.chained_assignment = None 


def permut2(x):
    xmax = x
    xout = x
    y = []

    while xout > 0:
        for xi in range(xout)[:-1]:
            y.append(int(xmax-xi)-1)
        xout-=1
    return y

def permut1(x):
    xmax = x
    xout = x
    y = []

    while xout > 0:
        for xi in range(xout)[:-1]:
            y.append(int(xmax-xout))
        xout-=1
    return y

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def func2(x, a, b, c): 
    return a*np.log(x + b)+c

def f1(x, a0, a1, a2, a3):
    return a0 + a1/x + a2/x**2 + a3/x**3






re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)
        
        
        
        
#version = 2
#print('Searching for file...')
path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\python scripts\\sizec ptmx\\'
filename = 'nSEC_atlas_v1_cleaned.xlsx'

df = pd.read_excel(path+filename)
print('File read.')
phmap = 'viridis'
CUTOFF = 1, 400

frag_rho = []
frag_pvals = []
frag_prot = []
frag_rhocon = []

#targets = [x for x in df.columns[2:]]
targets = ['CDK6']



for tg in targets:
    dfs = pd.DataFrame({'Model':df['Model'], 'Fraction':df['Fraction'], 'Target':df[tg]})
    dfs = dfs[dfs['Fraction']!=1].reset_index(drop=True)
    dfs = dfs[dfs['Fraction']!=24].reset_index(drop=True)

    dfscells = [x for x in dfs['Model'].unique()]
    dictresh = {'Fraction':dfs[dfs['Model']==dfscells[0]].reset_index(drop=True)['Fraction']}

    for cell in dfscells:
        dictresh[cell]=dfs[dfs['Model']==cell].reset_index(drop=True)['Target']

    dfr = pd.DataFrame(dictresh)

#    print(dfr.shape)
#    print(dfr.columns)

    if ([x for x in dfr.sum().values].count(0.0) < CUTOFF[0]+1) and (np.sum([x for x in dfr.sum().values]) > CUTOFF[1]):
        newlabels = [f"{cell} cellnorm" for cell in dfscells]
        dfrn = dfr.copy()
        for cell in dfscells:
            cellstat = sum(dfrn[cell])
            dfrn[f"{cell} cellnorm"] = dfrn[cell]/max(cellstat, 0.1)

        x = dfrn[newlabels].values
        pca = PCA(n_components=1)
        pCones = pca.fit_transform(x.T)
    
        dfpc = dfrn[newlabels]
        dfpc.loc[22, :] = [x[0] for x in pCones]
        dfpca = dfpc.sort_values(by=22, axis=1)
        dfpca=dfpca.drop(22)
        

        rhos = []
        pvals = []
        rhocon = []
        con = []
        lib = []
    
        for pmx, pmy in zip(permut1(len(dfpca.columns)), permut2(len(dfpca.columns))):
            rho, pval = stats.spearmanr(dfpca[dfpca.columns[pmx]],dfpca[dfpca.columns[pmy]])
            rhos.append(rho)
            rhocon.append(abs(rho))
            pvals.append(-np.log10(pval))
            con.append(dfpca.columns[pmx][:-9])
            lib.append(dfpca.columns[pmy][:-9])
            
        dfpmd = pd.DataFrame({'rho':rhos, 'pval':pvals, 'colL':con, 'colR':lib})
        dfpmd.to_excel(path+'spear_test.xlsx', index=False)
        rhomed = np.median(pd.Series(rhos).dropna())
        pmed = np.median(pd.Series(pvals).dropna())
        rhoconmed = np.median(pd.Series(rhocon).dropna())
        
        frag_rho.append(rhomed)
        frag_pvals.append(pmed)
        frag_prot.append(tg)
        frag_rhocon.append(rhoconmed)
        
        print('-', end='')
    
#fragdf = pd.DataFrame({'rho':frag_rho, 'rhocon':frag_rhocon,'pval':frag_pvals, 'prot':frag_prot})
#fragdf.to_excel(path+'fragdf_vmod.xlsx', index=False)
#    
    
#fragdf = pd.read_excel(path+'fragdf_vmod.xlsx')
#print('File read.')
#
#x = np.linspace(0, 1, 200)
#
#popt, pcov = curve_fit(func, fragdf['rhocon'], fragdf['pval'])
#perr = np.sqrt(np.diag(pcov))
#print(popt,'\n\n\n', pcov, '\n\n', perr, '\n\n', sum(perr))
#plt.scatter(fragdf['rhocon'], fragdf['pval'])
#plt.plot(x, func(x, *popt), c='green')
#plt.show()
#plt.clf()    
#
#fragdf['ph sep'] = fragdf['rhocon']-fragdf['rho']
#
#pmask = fragdf['pval']>-np.log10(0.05)
#smask = fragdf['ph sep']> fragdf['ph sep'].std()*5
#vmask = fragdf['pval']> fragdf['pval'].std()*5
#
#fragdf['hit'] = ''
#fragdf['rho lab'] = ''
#fragdf['sep lab'] = ''
#fragdf['lab'] = ''
#
##xf.loc[p_mask&u_mask3, f'{cat_code} Label {hit_std}std'] = xf[p_mask&u_mask3]['Gene Symbol']
##xf.loc[p_mask&d_mask3, f'{cat_code} Sig {hit_std}std'] = 'decreased expression'
#            
#fragdf.loc[pmask&smask, 'hit'] = 'sep'
#fragdf.loc[pmask&smask, 'sep lab'] = fragdf[pmask+smask]['prot']
#fragdf.loc[vmask, 'hit'] = 'rho'
#fragdf.loc[vmask, 'rho lab'] = fragdf[vmask]['prot']
#fragdf.loc[vmask, 'lab'] = fragdf[vmask]['prot']
#fragdf.loc[pmask&smask, 'lab'] = fragdf[pmask+smask]['prot']
#
#print(fragdf.shape)
#
#fragdf.to_excel(path+'frag_fin.xlsx', index=False)