# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:38:32 2024

@author: HUFFMP
"""
import pandas as pd
from proteomicRuler.ruler import Ruler, add_mw
# import sys
# sys.stdout = open('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\meth\\hmmm.txt', 'w')

switch='y'
while switch.lower()=='y':
    # path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\kat2global\\exports\\'
    # file = 'PFAS_PH_20231103_K22global_fx_prot.xlsx'

    pathfile = input('PD Output File Path: ')
    path = ('\\').join(pathfile.split('\\')[:-1])+'\\'
    file = pathfile.split('\\')[-1]

    df = pd.read_excel(pathfile)

    ploidy = float(input("Cell Ploidy (ex. 2): "))
    protconc = float(input("Cell Protein Concentration ug/uL (generally between 200-300): "))



    mw_id = 'Mass'
    mw_col = mw_id
    acc_id = "Protein IDs"

    df = add_mw(df, acc_id)
    print(df.columns)
    df = df[pd.notnull(df[mw_col])]
    df[mw_col] = df[mw_col].astype(float)

    ab_ids = [x for x in df.columns if str(x).count('Abundance:')>0]

    df = df.astype({mw_id: int})
    df = df.dropna(subset=ab_ids)

    print(f"{df.shape[0]} Proteins Quantified.")
    print("Calculating Copy Numbers.....")

    rdf = Ruler(df, ab_ids, mw_id, acc_id, ploidy, protconc)

    rdf.df.to_excel(path+file.split('.')[0]+'_Hruler.xlsx', index=False)
    switch = input('File saved. Analyze new file? (y/n): ')

print('exiting loop')
# sys.exit()
# sys.stdout.close()






























