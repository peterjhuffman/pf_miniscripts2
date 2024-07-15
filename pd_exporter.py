# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 09:36:00 2022

@author: HUFFMP
"""

import sys
import json
import os
import pandas as pd

json_loc = sys.argv[1]

with open(json_loc, encoding="utf-8") as g:
    lines = json.load(g)

infoloc1 = lines['Tables'][0]['DataFile']
infoloc2 = lines['Tables'][1]['DataFile']
infoloc3 = lines['Tables'][2]['DataFile']
infoloc4 = lines['Tables'][3]['DataFile']

jobname = lines['ResultFilePath'].split('.')[0].split('/')[-1]

df = pd.read_csv(infoloc1, sep='\t')
pdf = pd.read_csv(infoloc2, sep='\t')
sdf = pd.read_csv(infoloc3, sep='\t')
mdf = pd.read_csv(infoloc4, sep='\t')

newpath = f'\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\PD export\\{jobname}'
if not os.path.exists(newpath):
    os.makedirs(newpath)



df.to_excel(f'{newpath}\\{jobname}_prot.xlsx')
pdf.to_excel(f'{newpath}\\{jobname}_pep.xlsx')
sdf.to_excel(f'{newpath}\\{jobname}_psm.xlsx')
mdf.to_excel(f'{newpath}\\{jobname}_msms.xlsx')


