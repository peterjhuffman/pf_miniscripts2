# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:06:38 2024

@author: HUFFMP
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:34:40 2024

@author: HUFFMP
"""
import pandas as pd

filepath ='//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/results/37. INSM1 multidosing/exports/'
file = 'PFEC_PH_20240624_ISDglobal_fx_prot_pure_analysisA.xlsx'
kf = pd.read_csv(filepath+'uni itx.csv')

# knf=pd.DataFrame({'interactors':kf['interactors'].unique()})
# knf.to_csv(filepath+'uni itx.csv')


label = 'znf'

df = pd.read_excel(filepath+file)


known = [x for x in kf['interactors'].values]
print(known)

# xf.loc[bhc_mask&s0_mask, f'{cat_code} Label FDR'] = xf[bhc_mask&s0_mask]['Gene Symbol']
# xf.loc[bhc_mask&s0_mask&bhc_d_mask, f'{cat_code} Sig FDR'] = 'decreased expression'

kmask = df['Gene Symbol'].str.upper().isin(known)

df[f'{label} target'] = ''
df[f'{label} label'] = ''



df.loc[kmask, f'{label} substrate'] = 'target'
df.loc[kmask, f'{label} label'] = df[kmask]['Gene Symbol']

df.to_excel(filepath+file.split('.')[0]+'_itx.xlsx', index=False)




























