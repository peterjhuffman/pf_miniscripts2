# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:29:31 2024

@author: huffmp
"""


import pandas as pd
import gseapy as gp
import numpy as np
import matplotlib.pyplot as plt
import pathlib
from time import sleep, time
import sys
import os

def crash():
    print('\nExiting PHseus GSEA. \n\n')
    sleep(1.2)
    sys.exit()
    
    
def timer_quant(sec):
    """
    Splits a value in seconds into hours, minutes and seconds.
    Receives:
        sec : total seconds
    Returns:
        n_hr : number of hours
        n_min : number of minutes
        n_sec : number of seconds
    """
    nclock = round(sec, 1)
    n_hr = int(nclock/3600)
    n_min = int((nclock-(3600*n_hr))/60)
    n_sec = nclock - (3600*n_hr) - 60*n_min
    return (n_hr, n_min, n_sec)


def timed_print(message, timesplit):
    """
    Prints a message, plus the time since the last message sent.
    Receives:
        message: message to print
        timesplit : time of last message
    Returns:
        timesplit : current time
    """
    nclock = timer_quant(time()-timesplit)
    timesplit = time()
    print(message.ljust(55)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return timesplit


def enr_mod(df):
    df['-logp'] = -np.log10(df['Adjusted P-value'])
    del df['Old P-value']
    del df['Old Adjusted P-value']
    del df['Odds Ratio']

    df['hitcount'] = df['Genes'].str.count(';')
    df = df[df['hitcount']>0]
    del df['hitcount']
    return df


def write_enr_depict(df, path, note):
    df = df.sort_values(by='-logp', ascending=True)
    df = df[df['-logp']>(-np.log10(0.05))]
    
    set_tf = ['Transcription_Factor_PPIs']
    set_comp = ['CORUM']
    set_bio = ['KEGG_2021_Human', 'GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'WikiPathways_2019_Human']

    p_segments = [0, 2.6, 4, 7.2, 11, 18]

    colors_tf = {1:'#CDEDFE', 2:'#8ED7FD', 3:'#37B7FB', 4:'#0489D0', 5:'#02517B', 6:'#01273D'}
    colors_comp = {1:'#FECDED', 2:'#FD8ED7', 3:'#FB37B7', 4:'#D00489', 5:'#7B0251', 6:'#3D0127'}
    colors_bio = {1:'#EDFECD', 2:'#D7FD8E', 3:'#B7FB37', 4:'#89D004', 5:'#517B02', 6:'#273D01'}

    if df.shape[0]>0:
        plt.figure(figsize=(10,6))
        df_tf = df[df['Gene_set'].isin(set_tf)].reset_index(drop=True)
        df_comp = df[df['Gene_set'].isin(set_comp)].reset_index(drop=True)
        df_bio = df[df['Gene_set'].isin(set_bio)].reset_index(drop=True)

        if df_tf.shape[0]>0:
            df_tf = df_tf.head(10)
            cats = df_tf['Term']
            vals = df_tf['-logp']
            
            v_segments = [len([p for p in p_segments if p<=val]) for val in vals]
            colorplex = [colors_tf[v] for v in v_segments]
 
            plt.barh(cats, vals, color=colorplex)
            plt.xlabel('-log10(p)')
            plt.ylabel('transcription factor PPI')
            plt.title('Enriched Transcription Factors')
            plt.savefig(path+note+'_enrTF.png', bbox_inches="tight", dpi=300)
            plt.clf()
    
        if df_comp.shape[0]>0:
            df_comp = df_comp.head(10)
            cats = df_comp['Term']
            vals = df_comp['-logp']
            
            v_segments = [len([p for p in p_segments if p<=val]) for val in vals]
            colorplex = [colors_comp[v] for v in v_segments]
            
            plt.barh(cats, vals, color=colorplex)
            plt.xlabel('-log10(p)')
            plt.ylabel('protein complex')
            plt.title('Enriched Protein Complexes')
            plt.savefig(path+note+'_enrCOMP.png', bbox_inches="tight", dpi=300)
            plt.clf()
            
        if df_bio.shape[0]>0:
            df_bio = df_bio.head(10)
            cats = df_bio['Term']
            vals = df_bio['-logp']
            
            v_segments = [len([p for p in p_segments if p<=val]) for val in vals]
            colorplex = [colors_bio[v] for v in v_segments]
            
            plt.barh(cats, vals, color=colorplex)
            plt.xlabel('-log10(p)')
            plt.ylabel('biological pathway')
            plt.title('Pathway and Bionetwork Enrichment')
            plt.savefig(path+note+'_enrBIO.png', bbox_inches="tight", dpi=300)
            plt.clf()




#filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\30. k5ip\\exports\\'
#file = 'PFAS_PH_20240328_k5i_prot_pure_analysis2.xlsx'
#
#fp2 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\32. AHI insm1 rtsMS3 global\\exports\\'
#f2 = 'PFAS_PH_20240405_AHIglobal_rtsMS3_prot_analysis.xlsx'
#
#writepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\gseapy\\'
#
##df= pd.read_excel(filepath+file)
#df = pd.read_excel(fp2+f2)

filename = input('Protein Export Filename? : ')
if filename[0]=='f':
    filename = filename[5:]
fileSuffix = pathlib.Path(filename).suffix
foldersav = pathlib.Path(filename).stem
filepath = filename[:-(len(fileSuffix)+len(foldersav))]
print('Searching for and reading file.')
if fileSuffix == '.xlsx':
    df = pd.read_excel(filename)
elif fileSuffix == '.csv':
    df = pd.read_csv(filename)
else:
    print('File type not recognized.')
    crash()
print('File read successfully.')

clock = time()


condlist = []
complist = []

ab_count = 1
for ab in df.columns:
    if ab.count(' Sig ')>0:
        condlist = list(set(condlist + [ab.split(' Sig ')[0]]))
        complist = list(set(complist + [ab.split(' Sig ')[1]]))

newpath = f'{filepath}phGSEA_{foldersav}'
if not os.path.exists(newpath):
    os.makedirs(newpath)
print('Folders created. Analyzing enrichment: ')

for cond in condlist:
    condsav = cond.replace(' / ','_').replace(' ','')
    newpath = f'{filepath}phGSEA_{foldersav}\\{condsav}'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    for comp in complist:
        dhits = [x for x in df[df[f'{cond} Sig {comp}']=='decreased expression']['Gene Symbol'].dropna().values]
        uhits = [x for x in df[df[f'{cond} Sig {comp}']=='increased expression']['Gene Symbol'].dropna().values]



        if len(dhits)>0:
            enr = gp.enrichr(gene_list=dhits, # or "./tests/data/gene_list.txt",
                             gene_sets=['CORUM','Reactome_2022', 'KEGG_2021_Human',
                                        'GO_Biological_Process_2023', 'GO_Molecular_Function_2023',
                                        'WikiPathways_2019_Human', 'Transcription_Factor_PPIs'],
                             organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                             outdir=None, # don't write to disk
                            )
            
            ddf = enr.results
            ddf = enr_mod(ddf).reset_index(drop=True)
            write_enr_depict(ddf, newpath, f'\\downreg_{comp}_{condsav}')
            ddf.to_excel(newpath+f'\\downreg_{comp}_{condsav}_annotations.xlsx',index=False)
            print('-',end='')
#            sleep(60)
            
        if len(uhits)>0:
            enr = gp.enrichr(gene_list=uhits, # or "./tests/data/gene_list.txt",
                             gene_sets=['CORUM','Reactome_2022', 'KEGG_2021_Human',
                                        'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023',
                                        'WikiPathways_2019_Human', 'Transcription_Factor_PPIs'],
                             organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                             outdir=None, # don't write to disk
                            )
            udf = enr.results
            udf = enr_mod(udf).reset_index(drop=True)
            udf.to_excel(newpath+f'\\upreg_{comp}_{condsav}_annotations.xlsx',index=False)

            write_enr_depict(udf, newpath, f'\\upreg_{comp}_{condsav}')


#            sleep(60)
            print('-',end='')
print()
clock2 = timed_print('time elapsed: ', clock)





















