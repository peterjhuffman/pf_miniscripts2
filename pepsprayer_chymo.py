# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 13:23:47 2023

@author: HUFFMP
"""

"""
PEPsprayer v01

Model to estimate likelihood of peptide identification.

"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
from time import sleep as slp
from time import time as timer
import sys

pd.options.mode.chained_assignment = None

def pspray():
    pass


def fasta_interp(file, timesplit):
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Reading Fasta file.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    f = open(file, "r")
    lines = f.readlines()
    f.close()

    proteome = {}

    fasta = 'huffmp'
    code = 'admin'


    for line in lines:
        if line[0] == '>':
            proteome[code] = fasta
            code = line.split('|')[1].split(' ')[0]
            fasta = ''
        else:
            fasta += line[:-1]

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Fasta Interpreted.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return proteome, timesplit

def mol_define(seq):
    if len(seq)<900:
        return Chem.MolFromSmiles(seq)
    else:
        return 'x'


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

    if seq[-1] == 'y':
        smile = 'CC(=O)X'
        seq = seq[:-1]
    if seq[-1] == 'm':
        smile = 'CC(C)CN1CCCC1C(=O)NCCC(=O)NCCC(=O)X'
        seq = seq[:-1]
    else:
        smile = 'X'

    for char in seq:
        smile = smile[:-1]
        smile += aa_smiles[char]

    return smile


def translate_mods(seq, mods):
    seq = str(seq)
    if mods['C'] == 1:
        seq = seq.replace('C', 'c')
    if mods['O'] == 1:
        seq = seq.replace('M', 'o')
    if mods['TMTn'] == 1:
        seq += 'm'
    if mods['ACn'] == 1:
        seq += 'y'
    if mods['TMT']!=0:
        for tmt in mods['TMT']:
            seq = seq[:int(tmt)-1]+'t'+ seq[int(tmt):]
    if mods['AC']!=0:
        for ac in mods['AC']:
            seq = seq[:int(ac)-1]+'a'+ seq[int(ac):]
    return seq


def mod_interp(modstring):
    modstring = str(modstring)
    mods = {'TMTn':0, 'TMT':0, 'C':0, 'O':0, 'ACn':0, 'AC':0, 'dAm':0}
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
    if modstring.count('1xAcetyl [N-Term]')>0:
        mods['ACn'] = 1
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


def calc_target_descriptors(mol, val, tot):
    if mol == 'x':
        return {'fr_Nhpyrrole':0, 'fr_guanido':0, 'BCUT2D_CHGLO':0, 'MinEStateIndex':0, 'Chi4v':0, 'fr_ether':0, 'AvgIpc':0, 'fr_unbrch_alkane':0, 'BCUT2D_MWLOW':0, 'BCUT2D_MRHI':0, 'BCUT2D_LOGPHI':0, 'FpDensityMorgan1':0}
    if val%(tot/100)<1:
        print(f"{int(val//(tot/100))}%", end='     ')
    dx = {}

    dx['fr_Nhpyrrole'] = Descriptors.fr_Nhpyrrole(mol)
    dx['fr_guanido'] = Descriptors.fr_guanido(mol)
    dx['BCUT2D_CHGLO'] = Descriptors.BCUT2D_CHGLO(mol)
    dx['MinEStateIndex'] = Descriptors.MinEStateIndex(mol)
    dx['Chi4v'] = Descriptors.Chi4v(mol)
    dx['fr_ether'] = Descriptors.fr_ether(mol)
    dx['FpDensityMorgan1'] = Descriptors.FpDensityMorgan1(mol)
    dx['AvgIpc'] = Descriptors.AvgIpc(mol)
    dx['fr_unbrch_alkane'] = Descriptors.fr_unbrch_alkane(mol)
    dx['BCUT2D_MRHI'] = Descriptors.BCUT2D_MRHI(mol)
    dx['BCUT2D_MWLOW'] = Descriptors.BCUT2D_MWLOW(mol)
    dx['BCUT2D_LOGPHI'] = Descriptors.BCUT2D_LOGPHI(mol)
    return dx


def calc_target_descriptors_old(mol, val, tot):
    if val%(tot/100)<1:
        print(f"{int(val//(tot/100))}%", end='     ')
    dx = {}

    dx['EState_VSA4'] = Descriptors.EState_VSA4(mol)
    dx['fr_Ndealkylation2'] = Descriptors.fr_Ndealkylation2(mol)
    dx['MaxAbsEStateIndex'] = Descriptors.MaxAbsEStateIndex(mol)
    dx['AvgIpc'] = Descriptors.AvgIpc(mol)
    dx['NumSaturatedRings'] = Descriptors.NumSaturatedRings(mol)
    dx['VSA_EState9'] = Descriptors.VSA_EState9(mol)
    dx['FpDensityMorgan1'] = Descriptors.FpDensityMorgan1(mol)
    dx['fr_alkyl_halide'] = Descriptors.fr_alkyl_halide(mol)
    dx['FpDensityMorgan3'] = Descriptors.FpDensityMorgan3(mol)
    dx['BCUT2D_MRHI'] = Descriptors.BCUT2D_MRHI(mol)
    dx['BCUT2D_MWLOW'] = Descriptors.BCUT2D_MWLOW(mol)
    dx['MinEStateIndex'] = Descriptors.MinEStateIndex(mol)
    return dx


def calc_descriptors(df, timesplit):
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Peptides Written to File.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    chemlist = [Chem.MolFromSmiles(x) for x in df['SMILES string']]

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Peptides Rendered.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    tot = len(chemlist)
    descrs = [decsrs_logged(mol, val, tot) for val, mol in enumerate(chemlist)]
    dfd = pd.DataFrame(descrs)

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print()
    print('Molecular Descriptors calculated.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    return pd.concat([df, dfd], axis=1), timesplit

def timed_print(message, timesplit):
    """
    Prints a message, plus the time since the last message sent.
    Receives:
        message: message to print
        timesplit : time of last message
    Returns:
        timesplit : current time
    """
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print(message.ljust(55)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return timesplit


def calc_descriptors_2(df, timesplit):
    targets = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']
    timesplit = timed_print(f'{df.shape[0]} Peptides Written to File.', timesplit)

    chemlist = [mol_define(x) for x in df['SMILES string']]

    timesplit = timed_print(f'{df.shape[0]} Peptides Rendered.', timesplit)

    tot = len(chemlist)
    descrs = [calc_target_descriptors(mol, val, tot) for val, mol in enumerate(chemlist)]
    dfd = pd.DataFrame(descrs)

    print()
    timesplit = timed_print(f'{dfd.size} Molecular Descriptors calculated.', timesplit)

    return pd.concat([df, dfd], axis=1).dropna(subset=targets), timesplit


def seq2pep(df, timesplit):
    df = df.reset_index()
    modlist = [mod_interp(x) for x in df['Modifications']]

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Modifications Parsed.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    mod_seqlist = [translate_mods(seq, modlist[mod_id]) for mod_id, seq in enumerate(df['Sequence'])]

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Mods Converted to Structure.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    chemlist = [seq2mol(x) for x in mod_seqlist]

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Peptides Designed.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    df['SMILES string'] = pd.Series(chemlist)
    return df, timesplit


def binary_peps(df, contam_ids, proteome, timesplit):
    MODTYPE = 'Modifications'
#    MODTYPE = 'Modifications (all possible sites)'

    t_cols = [MODTYPE, 'Positions in Master Proteins', 'Annotated Sequence', 'Master Protein Accessions']

    df = df[t_cols]

    df['Sequence'] = [f.split('.')[1] for f in list(df['Annotated Sequence'])]
    del df['Annotated Sequence']

    df['Accession_fixed'] = [f.split(';')[0].split('-')[0].split(' ')[0][:6] for f in list(df['Master Protein Accessions'])]
    df['accession3'] = 0
    df.loc[df['Accession_fixed'].isin(contam_ids), 'accession3'] = 1
    df = df[df['accession3']==0].reset_index(drop=True)
    del df['Accession_fixed']
    del df['accession3']

    df['Accession'] = [x.split(';')[0] for x in df['Master Protein Accessions']]
    del df['Master Protein Accessions']

    prots = df['Accession'].unique()
    posmaster = df['Positions in Master Proteins'].unique()
    df['Acquired'] = 1

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('File trimmed and equalized.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    tmtlab = sum(df[MODTYPE].str.count('TMT').dropna()) > 0
    carbamido = sum(df[MODTYPE].str.count('Carbamido').dropna()) > 0

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Creating missed peptides.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    den = len(prots)
    nen = 1
    for prot in prots:
        print(f'Creating. prot {nen}/{den}: {prot}')
        nen+=1
        if prot in proteome.keys():
            fasta = proteome[prot]
        elif prot.split('-')[0] in proteome.keys():
            fasta = proteome[prot.split('-')[0]]
        else:
            continue
        bases = []
        i, j = 0, -1
        while i !=-1:
            i = min(fasta[j+1:].find('F'), fasta[j+1:].find('Y'), fasta[j+1:].find('W'))
            if i==-1:
                i = max(fasta[j+1:].find('F'), fasta[j+1:].find('Y'), fasta[j+1:].find('W'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()

        for val, base in enumerate(bases):
            pep1 = fasta[base+1:bases[min(len(bases)-1, val+1)]+1]
            pep2 = fasta[base+1:bases[min(len(bases)-1, val+2)]+1]

            label1, label2, lb1empty, lb2empty = '', '', True, True

            if carbamido:
                if pep1.count('C')>0:
                    c1 = [i for i, letter in enumerate(pep1) if letter == 'C']
                    label1 += f'{len(c1)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c1])}]'
                    lb1empty = False
                if pep2.count('C')>0:
                    c2 = [i for i, letter in enumerate(pep2) if letter == 'C']
                    label2 += f'{len(c2)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c2])}]'
                    lb2empty = False


            if tmtlab:
                if len([i for i, letter in enumerate(pep1) if ((letter == 'K') and(i<len(pep1)-1))])>0:
                    tmt1 = [i for i, letter in enumerate(pep1) if ((letter == 'K') and(i<len(pep1)-1))]
                    if not lb1empty:
                        label1 += '; '
                    label1 += f'{len(tmt1)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt1])}]'
                    lb1empty = False
                if len([i for i, letter in enumerate(pep2) if ((letter == 'K') and(i<len(pep2)-1))])>0:
                    tmt2 = [i for i, letter in enumerate(pep2) if ((letter == 'K') and(i<len(pep2)-1))]
                    if not lb2empty:
                        label2 += '; '
                    label2 += f'{len(tmt2)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt2])}]'
                    lb2empty = False
                if lb1empty:
                    label1 += '1xTMTpro [N-Term]'
                else:
                    label1 += '; 1xTMTpro [N-Term]'
                if lb2empty:
                    label2 += '1xTMTpro [N-Term]'
                else:
                    label2 += '; 1xTMTpro [N-Term]'

            pos1 = f"{prot} [{base+2}-{bases[min(len(bases)-1, val+1)]+1}]"
            pos2 = f"{prot} [{base+2}-{bases[min(len(bases)-1, val+2)]+1}]"

            pepline1 = [label1, pos1, pep1, prot, 0]
            pepline2 = [label2, pos2, pep2, prot, 0]

            addition_dex = [MODTYPE, 'Positions in Master Proteins', 'Sequence', 'Accession', 'Acquired']

            if pos1 not in posmaster:
                # df=df.append(pd.Series(pepline1, index=addition_dex), ignore_index=True)
                df = pd.concat([df, pd.DataFrame(pd.Series(pepline1, index=addition_dex)).transpose()], ignore_index=True)
            if pos2 not in posmaster:
                # df=df.append(pd.Series(pepline2, index=addition_dex), ignore_index=True)
                df = pd.concat([df, pd.DataFrame(pd.Series(pepline2, index=addition_dex)).transpose()], ignore_index=True)

    df = df[df['Sequence']!='']

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print(f'Missed peptides created. Writing to file'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return df, timesplit


def timer_quant(sec):
    nclock = round(sec, 1)
    n_hr = int(nclock/3600)
    n_min = int((nclock-(3600*n_hr))/60)
    n_sec = nclock - (3600*n_hr) - 60*n_min
    return (n_hr, n_min, n_sec)


def stop():
    print('Exiting pSpray.exe. Thanks for dropping by.')
    slp(2)
    sys.exit()


def simple_peps(df, contam_ids, proteome, timesplit):
    MODTYPE = 'Modifications'
#    MODTYPE = 'Modifications (all possible sites)'

    t_cols = [MODTYPE, 'Positions in Master Proteins', 'Annotated Sequence', 'Master Protein Accessions']#, 'RT [min] (by Search Engine): Sequest HT']

    df = df[t_cols]

    df['Sequence'] = [f.split('.')[1] for f in list(df['Annotated Sequence'])]
    del df['Annotated Sequence']

    df['Accession_fixed'] = [f.split(';')[0].split('-')[0].split(' ')[0][:6] for f in list(df['Master Protein Accessions'])]
    df['accession3'] = 0
    df.loc[df['Accession_fixed'].isin(contam_ids), 'accession3'] = 1
    df = df[df['accession3']==0].reset_index(drop=True)
    del df['Accession_fixed']
    del df['accession3']

    df['Accession'] = [x.split(';')[0] for x in df['Master Protein Accessions']]
    del df['Master Protein Accessions']

    prots = df['Accession'].unique()
    posmaster = df['Positions in Master Proteins'].unique()
#     df['Acquired'] = 1

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('File trimmed and equalized.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    tmtlab = sum(df[MODTYPE].str.count('TMT').dropna()) > 0
    carbamido = sum(df[MODTYPE].str.count('Carbamido').dropna()) > 0

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Creating missed peptides.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")

    den = len(prots)
    nen = 1
    for prot in prots:
        break
        print(f'Creating. prot {nen}/{den}: {prot}')
        nen+=1
        fasta = proteome[prot]
        bases = []
        i, j = 0, -1
        while i !=-1:
            i = min(fasta[j+1:].find('K'), fasta[j+1:].find('R'))
            if i==-1:
                i = max(fasta[j+1:].find('K'), fasta[j+1:].find('R'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()

        for val, base in enumerate(bases):
            pep1 = fasta[base+1:bases[min(len(bases)-1, val+1)]+1]
            pep2 = fasta[base+1:bases[min(len(bases)-1, val+2)]+1]

            label1, label2, lb1empty, lb2empty = '', '', True, True

            if carbamido:
                if pep1.count('C')>0:
                    c1 = [i for i, letter in enumerate(pep1) if letter == 'C']
                    label1 += f'{len(c1)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c1])}]'
                    lb1empty = False
                if pep2.count('C')>0:
                    c2 = [i for i, letter in enumerate(pep2) if letter == 'C']
                    label2 += f'{len(c2)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c2])}]'
                    lb2empty = False


            if tmtlab:
                if len([i for i, letter in enumerate(pep1) if ((letter == 'K') and(i<len(pep1)-1))])>0:
                    tmt1 = [i for i, letter in enumerate(pep1) if ((letter == 'K') and(i<len(pep1)-1))]
                    if not lb1empty:
                        label1 += '; '
                    label1 += f'{len(tmt1)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt1])}]'
                    lb1empty = False
                if len([i for i, letter in enumerate(pep2) if ((letter == 'K') and(i<len(pep2)-1))])>0:
                    tmt2 = [i for i, letter in enumerate(pep2) if ((letter == 'K') and(i<len(pep2)-1))]
                    if not lb2empty:
                        label2 += '; '
                    label2 += f'{len(tmt2)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt2])}]'
                    lb2empty = False
                if lb1empty:
                    label1 += '1xTMTpro [N-Term]'
                else:
                    label1 += '; 1xTMTpro [N-Term]'
                if lb2empty:
                    label2 += '1xTMTpro [N-Term]'
                else:
                    label2 += '; 1xTMTpro [N-Term]'

            pos1 = f"{prot} [{base+2}-{bases[min(len(bases)-1, val+1)]+1}]"
            pos2 = f"{prot} [{base+2}-{bases[min(len(bases)-1, val+2)]+1}]"

            pepline1 = [label1, pos1, pep1, prot, 0]
            pepline2 = [label2, pos2, pep2, prot, 0]

            addition_dex = [MODTYPE, 'Positions in Master Proteins', 'Sequence', 'Accession', 'Acquired']

            if pos1 not in posmaster:
                # df=df.append(pd.Series(pepline1, index=addition_dex), ignore_index=True)
                df = pd.concat([df, pd.DataFrame(pd.Series(pepline1, index=addition_dex)).transpose()], ignore_index=True)
            if pos2 not in posmaster:
                # df=df.append(pd.Series(pepline2, index=addition_dex), ignore_index=True)
                df = pd.concat([df, pd.DataFrame(pd.Series(pepline2, index=addition_dex)).transpose()], ignore_index=True)

    df = df[df['Sequence']!='']

    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print(f'Missed peptides created. Writing to file'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return df, timesplit


def master(filename):
    PATH = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\chymo TRPM8\\bESI training\\'
    clock = timer()
    
    timesplit = clock
    
    contamsfile = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\'+
                   'contaminants3.fasta')
    proteomefile = ('\\\\amer.pfizer.com\\pfizerfiles\\Research\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\'+
                    'homo sapiens.fasta')
    
    f = open(contamsfile, "r")
    lines = f.readlines()
    f.close()
    
    contam_toplines = []
    for line in lines:
        if line.count('>')>0:
            contam_toplines.append(line)
    
    contam_ids = []
    for line in contam_toplines:
        contam_ids.append(line[1:7].split(';')[0].split('-')[0].split(' ')[0])



    datasource = (PATH+filename)
    
    fdtype = datasource[-3:]
    if fdtype == 'csv':
        dframe = pd.read_csv(datasource)
    elif fdtype == 'lsx':
        dframe = pd.read_excel(datasource)
    else:
        print('filetype not recognized.')
        stop()
    
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('File read successfully.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    
    proteome_fasta, timesplit = fasta_interp(proteomefile, timesplit)
    
    bi_df, timesplit = binary_peps(dframe, contam_ids, proteome_fasta, timesplit)
#    bi_df, timesplit = simple_peps(dframe, contam_ids, proteome_fasta, timesplit)
    # bi_df.to_excel(PATH+filename[:-5]+'_binary.xlsx', index=False)
#    bi_df = pd.read_excel(PATH+filename[:-5]+'_binary.xlsx', index=False)
    
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('File read.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    
    #calc_descriptors(bi_df)
    
    mol_df, timesplit = seq2pep(bi_df, timesplit)
    #mol_df.to_excel(PATH+filename[:-5]+'_mol.xlsx', index=False)
    
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Peptides Created.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    
    desc_df, timesplit = calc_descriptors_2(mol_df, timesplit)
#    desc_df, timesplit = calc_descriptors(mol_df, timesplit)
    
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print('Descriptors Written to File.'.ljust(42)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    
    desc_df.to_excel(PATH+'desc\\'+filename[:-5]+'_desc.xlsx', index=False)
    
    nclock = timer_quant(timer()-clock)
    print(f"Runtime: {nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    print()
    print()

for x in [x for x in range(34,51)]+[x for x in range(67,101)]:
    print(f'Analyzing: PFAS_CHYpep_{str(x)}.xlsx')
    master(f'PFAS_CHYpep_{str(x)}.xlsx')

#master('ar5degRG_p02839_.csv')
