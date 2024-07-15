# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:05:15 2023

@author: HUFFMP
"""

import pandas as pd
import numpy as np
import sys
from time import sleep
import pathlib
from scipy.stats import ttest_ind as tti
import matplotlib.pyplot as plt

def permut2(x):
    xmax = x
    xout = x
    y = []

    while xout > 0:
        for xi in range(xout)[:-1]:
            y.append(str(xmax-xi))
        xout-=1
    return y

def power_law(x, a, b):
    return a * (x ** b)

def permut1(x):
    xmax = x
    xout = x
    y = []

    while xout > 0:
        for xi in range(xout)[:-1]:
            y.append(xmax-xout)
        xout-=1
    return y

def crash():
    print('\nExiting PHseus. \n\n')
    sleep(1.2)
    sys.exit()

print('Welcome to PHseus. The streamlined Proteomics Data Analysis Platform.')

rep = 6
replane = int(18/rep)
hit_ops = ['y', 'n', 'N', 'Y']

filename = input('Protein Export Filename? : ')
if filename[0]=='f':
    filename = filename[5:]
fileSuffix = pathlib.Path(filename).suffix
foldersav = pathlib.Path(filename).stem
print('Searching for and reading file.')
if fileSuffix == '.xlsx':
    df = pd.read_excel(filename)
elif fileSuffix == '.csv':
    df = pd.read_csv(filename)
else:
    print('File type not recognized.')
    crash()
print('File read successfully.')

rep = int(input('Replicates for this experiment? (1,2,3,6):'))
replane = int(18/rep)

print('Select Experiment Type:')
print('   1: Single Control')
print('   2: Double Control')
print('   3: All Comparisons')
print('   X: Exit')
exptype_in = input('This Experiment: ')

if exptype_in == '1':
    exptype = 1
elif exptype_in == '2':
    exptype = 2
elif exptype_in == '3':
    exptype = 3
elif exptype_in.upper() == 'X':
    crash()
else:
    print('Input not recognized.')
    crash()


ab_cols = [x for x in df.columns if str(x).count('Abundance:')>0]
if len(ab_cols) < 18:
    print('Insufficient abundance columns. Make sure you export your PD files with RAW abundance, not scaled.')
    crash()

else_cols = list(set(df.columns).difference(set(ab_cols)))
else_cols.sort()


categories = []
cat_timer = 1
while cat_timer <int(18/rep)+1:
    categories.append(input(f'   Lane Name {cat_timer}(X to exit): '))
    if categories[-1].upper()=='X':
        crash()
    cat_timer += 1
print()

if exptype<3:
    print(f'You have entered {exptype_in} control lane(s). Which lane is the first control?')
    ctrl = []
    ctrl.append(input('Lane #: '))
    if not (ctrl[0] in [str(x+1) for x in range(int(18/rep))]):
        print('Entry is not a number of one of the lanes.')
        crash()

if exptype == 2:
    print(f'You have entered {exptype_in} control lanes. Which lane is the second control?')
    ctrl.append(input('Second Lane #: '))
    if not (ctrl[1] in [str(x+1) for x in range(int(18/rep))]):
        print('Entry is not a number of one of the lanes.')
        crash()

if exptype==1:
    ctrl_loc = [ctrl[0] for x in categories]
    ctrl_loc[int(ctrl[0])-1] = '!'
    ctrl_loc_loc = [x for x in range(int(18/rep))]
elif exptype==2:
    ctrl_loc_loc = [x for x in range(int(18/rep))]
    print('Please assign each experimental lane to a corresponding control.')
    print(f'   1: {categories[int(ctrl[0])-1]}')
    print(f'   2: {categories[int(ctrl[1])-1]}')
    print('   X: Exit')
    print()
    ctrl_loc = []
    for i,j in enumerate(categories):
        if str(i+1) in ctrl:
            ctrl_loc.append('!')
        else:
            new_lane = input(f'Control Lane - {j}:')
            if new_lane in [str(x+1) for x in range(int(18/rep))]:
                ctrl_loc.append(new_lane)
            elif new_lane.upper() == 'X':
                crash()
            else:
                print('Input not recognized.')
                crash()
else:
    ctrl_loc = permut2(int(18/rep))
    ctrl_loc_loc = permut1(int(18/rep))
    ctrl = [str(x+1) for x in range(replane)]

print()
print(f'{df.shape[0]} Proteins Identified.')
df = df.dropna(subset=ab_cols)
for col in ab_cols:
    df = df[df[col]!=0.0]
print(f'{df.shape[0]} Proteins Quantified.')

psm_filter = input('Filter for >1 PSM? (Y/N):')
if not (psm_filter in hit_ops):
    print('Entry is not a valid option.')
    crash()
if psm_filter.upper()=='Y':
    df = df[df['# PSMs']>1]
    print(f'{df.shape[0]} Proteins After PSM Filter.')

uniq_filter = input('Filter for >1 Unique Peps? (Y/N):')
if not (uniq_filter in hit_ops):
    print('Entry is not a valid option.')
    crash()
if uniq_filter.upper()=='Y':
    df = df[df['# Unique Peptides']>1]
    print(f'{df.shape[0]} Proteins After UniqPep Filter.')


df_main = df[ab_cols].replace(0.0, 1.0)
df_aux = df[else_cols]


df_main_transformed = np.log2(df_main)

df_col_median = df_main_transformed.median()
xf = df_main_transformed-df_col_median

for col in else_cols:
    xf[col] = df_aux[col]

hit_std = 3
hit_d1 = input('Add Hit Detection - (p<0.05, diff>__*std)? (Y/N):')
if not (hit_d1 in hit_ops):
    print('Entry is not a valid option.')
    crash()
if hit_d1.upper()=='Y':
    hit_std = input('   Standard Deviation (Enter a value greater than 0)(Default 3):')
    if (float(hit_std)<0):
        print('Invalid stdev.')
        crash()
    hit_std = float(hit_std)


hit_fc = 2
hit_d3 = input('Add Hit Detection - (p<0.05, diff>__FC)? (Y/N):')
if not (hit_d3 in hit_ops):
    print('Entry is not a valid option.')
    crash()
if hit_d3.upper()=='Y':
    hit_fc = input('   Fold Change (Enter a value greater than 0)(Default 2):')
    if (float(hit_fc)<0):
        print('Invalid FC.')
        crash()
    hit_fc = float(hit_fc)


hit_d13 = input('Add Hit Detection - (p<0.05, diff>__FC AND diff>__*std)? (Y/N):')
if not (hit_d3 in hit_ops):
    print('Entry is not a valid option.')
    crash()
if hit_d13.upper()=='Y':
    hit_fc = input('   Fold Change (Enter a value greater than 0)(Default 2):')
    if (float(hit_fc)<0):
        print('Invalid FC.')
        crash()
    hit_fc = float(hit_fc)
    hit_std = input('   Standard Deviation (Enter a value greater than 0)(Default 3):')
    if (float(hit_std)<0):
        print('Invalid stdev.')
        crash()
    hit_std = float(hit_std)


hit_d4 = input('Add Hit Detection - (p<0.05, BH- FDR(custom))? (Y/N):')
if not (hit_d4 in hit_ops):
    print('Entry is not a valid option.')
    crash()
if hit_d4.upper()=='Y':
    fdr = input('   FDR (Enter a value between 0 and 1)(Default 0.05):')
    if (float(fdr)<0 or float(fdr)>1):
        print('Invalid FDR.')
        crash()
    fdr = float(fdr)
    s0 = input('   s0 (Enter a minimum LOG2fc for a hit if p=∞)(Default 0.2):')
    if (float(s0)<0):
        print('Invalid s0.')
        crash()
    s0 = float(s0)


hit_d5 = input('Add Hit Detection - (p<0.05, s0=0.33, FDR ex statuam)? (Y/N):')
if not (hit_d5 in hit_ops):
    print('Entry is not a valid option.')
    crash()



print('Calculating Fold Change...')
#print(ctrl_loc)
for i, loc in zip(ctrl_loc_loc, ctrl_loc):
    if loc != '!':
        cat_loc = categories[i]
        cat_code = f"{categories[i]} / {categories[int(ctrl[int(loc)-1])-1]}"
        col_bounds = (i*rep, (i*rep)+rep)
        df_cat = xf[ab_cols[col_bounds[0]:col_bounds[1]]]
        ctrl_bounds = ((int(ctrl[int(loc)-1])-1)*rep, (int(ctrl[int(loc)-1])-1)*rep+rep)
        df_ctrl = xf[ab_cols[ctrl_bounds[0]:ctrl_bounds[1]]]
        
        df_diff = df_cat.mean(axis=1)-df_ctrl.mean(axis=1)
        df_pval = -np.log10(tti(df_cat, df_ctrl, axis=1)[1])

        xf[f'{cat_code} -log(p)'] = df_pval
        xf[f'{cat_code} Difference'] = df_diff
        if exptype==3:
            xf[f'{cat_code} Difference REV'] = -df_diff

        p_cutoff = -np.log10(0.05)
        fc_cutoff = np.log2(hit_fc)

        p_mask = xf[f'{cat_code} -log(p)']>p_cutoff
        s_mask0 = xf[f'{cat_code} Difference'] > 0.33
        s_mask1 = xf[f'{cat_code} Difference'] < -0.33
        s_mask9 = xf[f'{cat_code} -log(p)']>power_law(np.abs(xf[f'{cat_code} Difference']), 2.42793373, -0.49068543)
        d_mask3 = xf[f'{cat_code} Difference'] < 0-xf[f'{cat_code} Difference'].std()*hit_std
        u_mask3 = xf[f'{cat_code} Difference'] > xf[f'{cat_code} Difference'].std()*hit_std
        d_mask0 = xf[f'{cat_code} Difference'] < -1*fc_cutoff
        u_mask0 = xf[f'{cat_code} Difference'] > fc_cutoff
        bhc_u_mask = xf[f'{cat_code} Difference'] > 0
        bhc_d_mask = xf[f'{cat_code} Difference'] < 0

        if hit_d1.upper() == 'Y':
            xf[f'{cat_code} Sig {hit_std}std'] = ''
            xf[f'{cat_code} Label {hit_std}std'] = ''
            xf.loc[p_mask&d_mask3, f'{cat_code} Label {hit_std}std'] = xf[p_mask&d_mask3]['Gene Symbol']
            xf.loc[p_mask&u_mask3, f'{cat_code} Label {hit_std}std'] = xf[p_mask&u_mask3]['Gene Symbol']
            xf.loc[p_mask&d_mask3, f'{cat_code} Sig {hit_std}std'] = 'decreased expression'
            xf.loc[p_mask&u_mask3, f'{cat_code} Sig {hit_std}std'] =  'increased expression'

        if hit_d3.upper() == 'Y':
            xf[f'{cat_code} Sig {hit_fc}FC'] = ''
            xf[f'{cat_code} Label {hit_fc}FC'] = ''
            xf.loc[p_mask&d_mask0, f'{cat_code} Label {hit_fc}FC'] = xf[p_mask&d_mask0]['Gene Symbol']
            xf.loc[p_mask&u_mask0, f'{cat_code} Label {hit_fc}FC'] = xf[p_mask&u_mask0]['Gene Symbol']
            xf.loc[p_mask&d_mask0, f'{cat_code} Sig {hit_fc}FC'] = 'decreased expression'
            xf.loc[p_mask&u_mask0, f'{cat_code} Sig {hit_fc}FC'] =  'increased expression'

        if hit_d13.upper() == 'Y':
            xf[f'{cat_code} Sig {hit_fc}FC {hit_std}std'] = ''
            xf[f'{cat_code} Label {hit_fc}FC {hit_std}std'] = ''
            xf.loc[p_mask&d_mask0&d_mask3, f'{cat_code} Label {hit_fc}FC {hit_std}std'] = xf[p_mask&d_mask0]['Gene Symbol']
            xf.loc[p_mask&u_mask0&u_mask3, f'{cat_code} Label {hit_fc}FC {hit_std}std'] = xf[p_mask&u_mask0]['Gene Symbol']
            xf.loc[p_mask&d_mask0&d_mask3, f'{cat_code} Sig {hit_fc}FC {hit_std}std'] = 'decreased expression'
            xf.loc[p_mask&u_mask0&u_mask3, f'{cat_code} Sig {hit_fc}FC {hit_std}std'] =  'increased expression'

        if hit_d4.upper() == 'Y':
            df_bhc = pd.DataFrame({'p':df_pval, 'dif':abs(df_diff),'index':[x for x in range(len(df_pval))]})
            df_bhc = df_bhc.sort_values(by='p', ascending=False)
            df_bhc['rank'] = [x+1 for x in range(len(df_pval))]
            df_bhc['bhc'] = [fdr*(x/len(df_pval)) for x in df_bhc['rank']]
            df_bhc['p adj'] = 10**(-df_bhc['p'])
            
            bhc_mask0 = df_bhc['rank']<=df_bhc[df_bhc['p adj']<df_bhc['bhc']]['rank'].max()
            
            df_bhc['sig'] = ''
            df_bhc.loc[bhc_mask0, 'sig'] = '+'

            df_bhc = df_bhc.sort_values(by='index')

            bhc_mask = df_bhc['sig']=='+'
            s0_mask = abs(xf[f'{cat_code} Difference']) > s0

            xf[f'{cat_code} Sig FDR'] = ''
            xf[f'{cat_code} Label FDR'] = ''
            xf.loc[bhc_mask&s0_mask, f'{cat_code} Label FDR'] = xf[bhc_mask&s0_mask]['Gene Symbol']
            xf.loc[bhc_mask&s0_mask, f'{cat_code} Label FDR'] = xf[bhc_mask&s0_mask]['Gene Symbol']
            xf.loc[bhc_mask&s0_mask&bhc_d_mask, f'{cat_code} Sig FDR'] = 'decreased expression'
            xf.loc[bhc_mask&s0_mask&bhc_u_mask, f'{cat_code} Sig FDR'] =  'increased expression'

        if hit_d5.upper() == 'Y':
            xf[f'{cat_code} Sig statuam'] = ''
            xf[f'{cat_code} Label statuam'] = ''
            xf.loc[p_mask&s_mask0&s_mask9, f'{cat_code} Label statuam'] = xf[p_mask&s_mask0&s_mask9]['Gene Symbol']
            xf.loc[p_mask&s_mask1&s_mask9, f'{cat_code} Label statuam'] = xf[p_mask&s_mask1&s_mask9]['Gene Symbol']
            xf.loc[p_mask&s_mask1&s_mask9, f'{cat_code} Sig statuam'] = 'decreased expression'
            xf.loc[p_mask&s_mask0&s_mask9, f'{cat_code} Sig statuam'] =  'increased expression'

# lab = '2 / 1'
# xfa = xf[xf[lab+' Sig FDR']=='']
# xfb = xf[xf[lab+' Sig FDR']!='']
# plt.scatter(xfa[lab+' Difference'], xfa[lab+' -log(p)'])
# plt.scatter(xfb[lab+' Difference'], xfb[lab+' -log(p)'])
# plt.show()

print(f"Saving to {filename[:-len(fileSuffix)]}_analysis.xlsx")
xf.to_excel(f"{filename[:-len(fileSuffix)]}_analysis.xlsx", index=False)
# df_bhc.to_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\ubr5global\\bruh.xlsx')
input('PHseus Complete. Feel free to close window.')