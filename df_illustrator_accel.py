# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:37:13 2024

@author: HUFFMP
"""

import pandas as pd
import numpy as np
import sys
from time import sleep as slp



def stop():
    print('Exiting pd_illustrator.py. Thanks for dropping by.')
    slp(2)
    sys.exit()

datasource = input('File Location (x to exit):')
if datasource.lower()=='x':
    stop()
fdtype = datasource[-3:]
if fdtype == 'csv':
    df = pd.read_csv(datasource)
elif fdtype == 'lsx':
    df = pd.read_excel(datasource)
else:
    print('filetype not recognized.')
    stop()
print('File read successfully.')

ab_count = 1
for ab in df.columns:
    print(f"{ab_count}: {ab} {' '*(40-len(f'{ab_count}: {ab}'))}  Example:{df.iloc[1][ab_count-1]}")
    ab_count+=1
print()

ab_select = input('Enter the index of the SWISSPROT Accesion ID column for analysis.\n'+
                  'Example: if column 4 is desired, input should be: 4\nIndex (x to exit): ')
if ab_select.lower()=='x':
    stop()

sapien


























