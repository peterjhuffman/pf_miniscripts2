# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 09:36:00 2022

@author: HUFFMP
"""

import sys
import json
import os
import pandas as pd
import datetime

BSAcode = 'A0A140T897'
json_loc = sys.argv[1]
ms_machine = sys.argv[2]
ms_standard = sys.argv[3]

with open(json_loc, encoding="utf-8") as g:
    lines = json.load(g)

infoloc1 = lines['Tables'][0]['DataFile']
infoloc2 = lines['Tables'][1]['DataFile']
jobname = lines['ResultFilePath'].split('.')[0].split('/')[-1]

cols = ['initial', 'date', 'standard', 'gradient', 'injection,ng/fm', 'protein', 'peptides', 'PSMs', 'MS2', 'sequence coverage for BSA', 'filename']

df1 = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\ms_standardlog.xlsx', 'Lumos')
df2 = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\ms_standardlog.xlsx', 'Eclipse')
df3 = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\ms_standardlog.xlsx', 'Ascend')

df = pd.read_csv(infoloc1, sep='\t')
sdf = pd.read_csv(infoloc2, sep='\t')
user = os.getlogin()
date = str(datetime.date.today())

t_col_coverage = ''
for b in df.columns:
    if b.lower().count('coverage')>0:
        t_col_coverage = b


if ms_standard.lower() == 'bsa' and t_col_coverage:
    bsa_cov = df[df['Accession']==BSAcode].reset_index().loc[0][t_col_coverage]
    gradient = '0.5h'
    injection = '60fm'
    additive = pd.Series([user, date, 'BSA', gradient, injection, '', '', '', '', bsa_cov, jobname], index=cols)
    t_col_coverage = True

elif ms_standard.lower() == 'hela':
    gradient = '2h'
    injection = '200ng'
    prot=sdf[sdf['Name']=='Protein Groups - # Peptides'].reset_index().loc[0]['Count']
    pep=sdf[sdf['Name']=='Peptide Groups - # Proteins'].reset_index().loc[0]['Count']
    psm=sdf[sdf['Name']=='PSMs - # Proteins'].reset_index().loc[0]['Count']
    ms2=sdf[sdf['Name']=='MS/MS Spectrum Info - # Precursors'].reset_index().loc[0]['Count']
    additive = pd.Series([user, date, 'Hela', gradient, injection, prot, pep, psm, ms2, '', jobname], index=cols)
    t_col_coverage = True

if t_col_coverage:
    if ms_machine.lower() == 'lumos':
        df1 = df1.append(additive, ignore_index=True)
    elif ms_machine.lower() == 'eclipse':
        df2 = df2.append(additive, ignore_index=True)
    elif ms_machine.lower() == 'ascend':
        df3 = df3.append(additive, ignore_index=True)
    
    with pd.ExcelWriter('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\ms_standardlog.xlsx') as writer:
        df1.to_excel(writer, sheet_name='Lumos', index=False)
        df2.to_excel(writer, sheet_name='Eclipse', index=False)
        df3.to_excel(writer, sheet_name='Ascend', index=False)
