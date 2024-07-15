# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:32:48 2023

@author: HUFFMP
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
import sklearn
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle

pd.options.mode.chained_assignment = None  # default='warn'

def quantize(x):
    return int(round(x, 0))


def ifprint(y, m):
    if y:
        print(m)

def translate_mods(seq, mods):
    yn = False
    seq = str(seq)
    ifprint(yn, seq)
    ifprint(yn, seq)
    if mods['C'] == 1:
        seq = seq.replace('C', 'c')
    ifprint(yn, seq)
    if mods['O'] == 1:
        seq = seq.replace('M', 'o')
    ifprint(yn, seq)
    if mods['TMTn'] == 1:
        seq += 'm'
    ifprint(yn, seq)
    if mods['ACn'] == 1:
        seq += 'y'
    ifprint(yn, seq)
    if mods['TMT']!=0:
        for tmt in mods['TMT']:
            seq = seq[:int(tmt)-1]+'t'+ seq[int(tmt):]
    ifprint(yn, seq)
    if mods['AC']!=0:
        for ac in mods['AC']:
            seq = seq[:int(ac)-1]+'a'+ seq[int(ac):]
    if mods['Ml'] == 1:
        seq = seq[1:]
    ifprint(yn, seq)
    return seq


def mod_interp(modstring):
    modstring = str(modstring)
    mods = {'TMTn':0, 'TMT':0, 'C':0, 'O':0, 'ACn':0, 'AC':0, 'dAm':0, 'Ml':0}
    if (modstring.count('1xTMTpro [N-Term]')+modstring.count('1xTMTpro [N-Term]'))>0:
        mods['TMTn'] = 1
    if modstring.count('xTMTpro [K')>0:
        mods['TMT'] = modstring.split('xTMTpro [K')[1].split(']')[0].split('; K')
        if mods['TMT'].count('')>0:
            mods['TMT'].remove('')
        if mods['TMT']==[]:
            mods['TMT']=0
    if modstring.count('xCarbamidomethyl')>0:
        mods['C'] = 1
    if modstring.count('xOxidation')>0:
        mods['O'] = 1
    if modstring.count('Acetyl [N-Term]')>0:
        mods['ACn'] = 1
    if modstring.count('xMet-loss')>0:
        mods['Ml'] = 1
    if modstring.count('xAcetyl [K')>0:
        mods['AC'] = modstring.split('xAcetyl [K')[1].split(']')[0].split('; K')
    if modstring.count('xDeamidated ')>0:
        mods['dAm'] = [i for i in modstring.split('xDeamidated [')[1].split(']')[0].split('; ') if i.count('/')==0]
    if mods['dAm'] == []:
        mods['dAm'] = 0
#    print(mods)
    return mods

def decsrs_logged(mol, val, tot):
    if val%(tot/100)<1:
        print(f"{int(val//(tot/100))}%", end='     ')
    return Descriptors.CalcMolDescriptors(mol)


def seq2pep(df):
    modlist = [mod_interp(x) for x in df['Modifications']]

    mod_seqlist = [translate_mods(seq, modlist[mod_id]) for mod_id, seq in enumerate(df['Sequence'])]

    chemlist = [seq2mol(x) for x in mod_seqlist]

    df['modseq'] = mod_seqlist
    df['SMILES string'] = pd.Series(chemlist)

    chemlist = [Chem.MolFromSmiles(x) for x in df['SMILES string']]
    df['MW'] = [Chem.Descriptors.MolWt(x) for x in chemlist]



#    descrs = [Descriptors.CalcMolDescriptors(mol) for mol in chemlist]
#    tot = len(chemlist)
#    descrs = [decsrs_logged(mol, val, tot) for val, mol in enumerate(chemlist)]

#    ddf = pd.DataFrame(descrs)
#    for col in ddf.columns:
#        df[col] = ddf[col]

    df['z (pred)'] = df['MW']/df['m/z [Da] (by Search Engine): Sequest HT']


    return df


def seq2mol(seq):
    aa_smiles = {'A': 'NC(C)C(=O)O',
                 'C': 'NC(CS)C(=O)O', 
                 'D': 'NC(CC(=O)O)C(=O)O',
                 'E': 'NC(CCC(=O)O)C(=O)O',
                 'F': 'NC(Cc1ccccc1)C(=O)O',
                 'G': 'NCC(=O)O',
                 'H': 'NC(Cc1c[nH]cn1)C(=O)O',
                 'I': 'NC(C(CC)C)C(=O)O',
                 'K': 'NC(CCCCN)C(=O)O',
                 'L': 'NC(CC(C)(C))C(=O)O',
                 'M': 'NC(CCSC)C(=O)O',
                 'N': 'NC(CC(=O)(N))C(=O)O',
                 'P': 'N1CCCC1C(=O)O',
                 'Q': 'NC(CCC(=O)(N))C(=O)O',
                 'R': 'NC(CCCNC(=N)(N))C(=O)O',
                 'S': 'NC(CO)C(=O)O',
                 'T': 'NC(C(O)C)C(=O)O',
                 'U': 'NC(C[SeH])C(=O)O',
                 'V': 'NC(C(C)C)C(=O)O',
                 'W': 'NC(Cc1c[nH]c2ccccc12)C(=O)O',
                 'Y': 'NC(Cc1ccc(O)cc1)C(=O)O',
                 't': 'NC(CCCCNC(=O)CCNC(=O)CCNC(=O)C1CCCN1(CC(C)C))C(=O)O',
                 'a': 'NC(CCCCNC(=O)C)C(=O)O', 
                 'o': 'NC(CCS(=O)C)C(=O)O', 
                 'c': 'NC(CSCC(=O)N)C(=O)O'}

    sw = False
    smile = 'X'
    if seq[-1] == 'y':
        smile = 'CC(=O)X'
        seq = seq[:-1]
    elif seq[-1] == 'm':
        smile = 'CC(C)CN1CCCC1C(=O)NCCC(=O)NCCC(=O)X'
        seq = seq[:-1]
    if seq[-1] == 'm':
        smile = 'CC(C)CN1CCCC1C(=O)NCCC(=O)NCCC(=O)N(C(=O)C)X'
        sw = True
        seq = seq[:-1]

    for char in seq:
        smile = smile[:-1]
        if sw:
            smile += aa_smiles[char][1:]
            sw = False
        else:
            smile += aa_smiles[char]

    return smile

aas = ['A', 'D', 'E','F','G','H','I', 'L','M','N','P','Q','R','S','T','V','W','Y', 'o', 't', 'c', 'y']


filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\'
df = pd.read_excel(filepath+'PFAS_PH_20231103_K22global_fx_mzpred.xlsx')

cat_col = ['Annotated Sequence','Modifications',
           'm/z [Da] (by Search Engine): Sequest HT']

denom = len(df[cat_col[0]])
xf = df[cat_col]
#print(xf.head())

seq = [str(x).split('.')[1] for x in xf[cat_col[0]]]
xf['Sequence'] = seq


mdf = seq2pep(xf)
mdf.to_excel(filepath+'_mdf.xlsx',index=False)

#aas = ['A', 'C', 'D', 'E','F','G','H','I', 'K', 'L','M','N','P','Q','R','S','T','U','V','W','Y', 'o', 't', 'a', 'c', 'y', 'm']
#aas = ['H', 'L','P','Q','R', 't']
#aas = ['K', 'R','H','o', 't', 'a', 'c', 'y', 'm']

#dkeys= Descriptors.CalcMolDescriptors(Chem.MolFromSmiles('C')).keys()
dkeys = ['MW']

AAc = pd.DataFrame()
for aa in aas:
    AAc[aa] = mdf['modseq'].str.count(aa)
    mdf[aa] = mdf['modseq'].str.count(aa)
#for desc in dkeys:
#    AAc[desc] = mdf[desc]

X = AAc
y = mdf['z (pred)']

#regr =sklearn.linear_model.LinearRegression()
#regr.fit(X, y)
#
#with open(filepath+'linreg.pkl', 'wb') as fid:
#    pickle.dump(regr, fid)

with open(filepath+'linreg.pkl', 'rb') as fid:
    regr = pickle.load(fid)
    regr.fit(X, y)

mdf['z (mod) lin'] = [regr.predict(np.array([[x.count(aa) for aa in aas]]))[0] for val,x in enumerate(mdf['modseq'])]
# +[y for y in mdf.iloc[val][dkeys].values]


score = sum([(x-y)**2 for y,x in zip(mdf['z (mod) lin'], mdf['z (pred)'])])
print(round(100*(score/denom), 1))

#mdf.to_excel(filepath+'_theomol.xlsx',index=False)

mdf['z (mod) lin quant'] = [quantize(x) for x in mdf['z (mod) lin']]
score = sum([(x-y)**2 for y,x in zip(mdf['z (mod) lin quant'], mdf['z (pred)'])])
print(round(100*(score/denom), 1))


mdf['z (pred) quant'] = [round(x,0)for x in mdf['z (pred)']]
mdf['check'] = [int(x) for x in mdf['z (pred) quant']==mdf['z (mod) lin quant']]

#mdf.to_excel(filepath+'_theomol.xlsx',index=False)
for val,x in enumerate(AAc.columns):
    print(AAc.columns[val], regr.coef_[val])
print(regr.intercept_)

print(round(100*mdf['check'].sum()/len(mdf['check']), 1))

#for i in range(len(aas)):
#    n_aas = aas.copy()
#    n_aas.pop(i)
#    print(n_aas)
#    print(aas[i])
#
#    AAc = pd.DataFrame()
#    for aa in n_aas:
#        AAc[aa] = mdf['modseq'].str.count(aa)
#    
#    X = AAc
#    y = mdf['z (pred)']
#    
#    
#    
#    regr =sklearn.linear_model.LinearRegression()
#    regr.fit(X, y)
#    
#    mdf['z (mod) lin'] = [regr.predict(np.array([[x.count(aa) for aa in n_aas]]))[0] for x in mdf['modseq']]
#    
#    score = sum([(x-y)**2 for y,x in zip(mdf['z (mod) lin'], mdf['z (pred)'])])
#    print(score)
#    
#    #mdf.to_excel(filepath+'_theomol.xlsx',index=False)
#    
#    mdf['z (mod) lin quant'] = [quantize(x) for x in mdf['z (mod) lin']]
#    score = sum([(x-y)**2 for y,x in zip(mdf['z (mod) lin quant'], mdf['z (pred)'])])
#    print(score)
#    print()


