# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:32:48 2023

@author: HUFFMP
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle
from time import sleep
import sys
import pathlib

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


def crash():
    print('\nExiting MZpred. \n\n')
    sleep(1.2)
    sys.exit()


aas = ['A', 'D', 'E','F','G','H','I', 'L','M','N','P','Q','R','S','T','V','W','Y', 'o', 't', 'c', 'y']


filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\'

print('Select Submission Type:')
print('   1: Single Sequence')
print('   2: Batch File Upload')
print('   X: Exit')
exptype_in = input('This Upload: ')
if exptype_in == '1':
    exptype = 1
elif exptype_in == '2':
    exptype = 2
elif exptype_in.upper() == 'X':
    crash()
else:
    print('Input not recognized.')
    crash()


if exptype == 1:
    filename = input('Peptide Sequence? :')
    

if exptype == 2:
    filename = input('Protein Export Filename? : ')
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

    ab_count = 1
    for ab in df.columns:
        print(f"{ab_count}: {ab} {' '*(70-len(f'{ab_count}: {ab}'))}  Example:{df.iloc[1][ab_count-1]}")
        ab_count+=1
    print()
    
    cat_col = input('Enter the index of the SEQUENCE or ANNOTATED SEQUENCE column for analysis.\n'+
                      'Example: if column 4 is desired, input should be: 4\nIndex (x to exit): ')
    if cat_col.lower()=='x':
        crash()

    cat_col = df.columns[int(cat_col)-1]


    xf = df[[cat_col]]
    #print(xf.head())


    if len(str(xf[[cat_col]].values[0]).split('.'))>1:
        seq = [str(x).split('.')[0] for x in xf[[cat_col]]]
    else:
        seq = [str(x).split('.')[1] for x in xf[[cat_col]]]


    xf['Sequence'] = seq
    mdf = seq2pep(xf)
    
    aas = ['A', 'C', 'D', 'E','F','G','H','I', 'K', 'L','M','N','P','Q','R','S','T','U','V','W','Y', 'o', 't', 'a', 'c', 'y', 'm']
    #aas = ['H', 'L','P','Q','R', 't']
    #aas = ['K', 'R','H','o', 't', 'a', 'c', 'y', 'm']
    
    
    AAc = pd.DataFrame()
    for aa in aas:
        AAc[aa] = mdf['modseq'].str.count(aa)
        mdf[aa] = mdf['modseq'].str.count(aa)
    #for desc in dkeys:
    #    AAc[desc] = mdf[desc]
    
    X = AAc
    
    #regr =sklearn.linear_model.LinearRegression()
    #regr.fit(X, y)
    #
    #with open(filepath+'linreg.pkl', 'wb') as fid:
    #    pickle.dump(regr, fid)
    
    with open(filepath+'linreg.pkl', 'rb') as fid:
        regr = pickle.load(fid)

    mdf['z (mod) lin'] = [regr.predict(np.array([[x.count(aa) for aa in aas]]))[0] for x in mdf['modseq']]
    mdf['z (mod) lin quant'] = [quantize(x) for x in mdf['z (mod) lin']]



    print(mdf.columns)



    print(foldersav)
#    mdf.to_excel(foldersav)


