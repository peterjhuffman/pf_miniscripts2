# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 12:19:13 2023

@author: HUFFMP
"""

import numpy as np
import pandas as pd

fp = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\pd concat\\'


fils = ['prot', 'pep', 'psms']

#t1 = pd.read_excel(fp+f'concat_10_{fils[0]}.xlsx')
#t2 = pd.read_excel(fp+f'concat_10_{fils[1]}.xlsx')
#t3 = pd.read_excel(fp+f'concat_10_{fils[2]}.xlsx')
#
#d1 = pd.read_excel(fp+f'concat_12_{fils[0]}.xlsx')
#d2 = pd.read_excel(fp+f'concat_12_{fils[1]}.xlsx')
#d3 = pd.read_excel(fp+f'concat_12_{fils[2]}.xlsx')
#
#x1 = pd.read_excel(fp+f'concat_2fx_{fils[0]}.xlsx')
#x2 = pd.read_excel(fp+f'concat_2fx_{fils[1]}.xlsx')
#x3 = pd.read_excel(fp+f'concat_2fx_{fils[2]}.xlsx')

y3 = t3.append(d3).reset_index(drop=True)


for col in t2.columns:
    if col.count('Found in Sample:')>1:
        del t2[col]
for col in d2.columns:
    if col.count('Found in Sample:')>1:
        del d2[col]

int2 = t2.append(d2).reset_index(drop=True)
int2['pepid temp'] = [str(x)+str(y) for x,y in zip(int2['Positions in Master Proteins'],int2['Modifications'])]

y2 = pd.DataFrame(columns=t2.columns)

int2['pepid count'] = [len(int2[int2['pepid temp']==x]['pepid temp']) for x in int2['pepid temp']]

int2.to_csv(fp+'int2.csv',index=False)


#for pepid in int2['pepid temp']:
#    sl = int2[int2['pepid temp']==pepid]
#    if sl.shape[0]==1:
#        y2 = y2.append(sl)
#    else:
#        

















































