# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:28:00 2023

@author: HUFFMP
"""

from sklearn.model_selection import train_test_split
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import confusion_matrix
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.decomposition import PCA
from time import time as timer
from random import sample
import pickle
import sys
from time import sleep as slp
import requests
import re
from requests.adapters import HTTPAdapter, Retry

pd.options.mode.chained_assignment = None



re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def timer_quant(sec):
    nclock = round(sec, 1)
    n_hr = int(nclock/3600)
    n_min = int((nclock-(3600*n_hr))/60)
    n_sec = nclock - (3600*n_hr) - 60*n_min
    return (n_hr, n_min, n_sec)


def timed_print(message, timesplit):
    nclock = timer_quant(timer()-timesplit)
    timesplit = timer()
    print(message.ljust(55)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return timesplit


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def stop():
    print('Exiting protein_ESI.py. Thanks for dropping by.')
    slp(2)
    sys.exit()


def read_model(filename):
    with open(filename, 'rb') as fid:
        clf_R = pickle.load(fid)
    return clf_R


def uniprot_call(acc):
    url = f'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+(organism_id:9606)+AND+accession:{str(acc)}&format=fasta'
    for batch, total in get_batch(url):
        d1 = {}
        d1['header'] = batch.text.splitlines()[0]
        d1['fasta'] = ''.join(batch.text.splitlines()[1:])
    return d1


def silico_digest(prot, tmtlab, alkyl, enzym, mcMAX):
    carbamido = ((alkyl=='IAA') or (alkyl=='CAA'))
    mcMAXt = mcMAX

    peps = pd.DataFrame(columns=['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites'])

    if enzym['Trypsin']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('K'), prot[j+1:].find('R'))
            if i==-1:
                i = max(prot[j+1:].find('K'), prot[j+1:].find('R'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)
    
        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'Trypsin', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)


    if enzym['GluC']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('D'), prot[j+1:].find('E'))
            if i==-1:
                i = max(prot[j+1:].find('D'), prot[j+1:].find('E'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)

        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'GluC', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)



    if enzym['ArgC']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('R'))
            if i==-1:
                i = max(prot[j+1:].find('R'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)
    
        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'ArgC', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)



    if enzym['ProA']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('P'), prot[j+1:].find('A'))
            if i==-1:
                i = max(prot[j+1:].find('P'), prot[j+1:].find('A'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)
    
        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'ProA', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)


    if enzym['Chymotrypsin']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('Y'), prot[j+1:].find('W'), prot[j+1:].find('F'))
            if i==-1:
                i = max(prot[j+1:].find('Y'), prot[j+1:].find('W'), prot[j+1:].find('F'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)
    
        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'Chymotrypsin', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)


    if enzym['LysC']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = min(prot[j+1:].find('K'))
            if i==-1:
                i = max(prot[j+1:].find('K'))
            bases.append(i+j+1)
            j += i+1
        bases.pop()
        bases.append(70000)
    
        for val, base in enumerate(bases):
            peplist, labels, lbempty = [], [], []
            mcMAXt = mcMAX
            while mcMAXt > 0:
                peplist.append(prot[base+1:bases[min(len(bases)-1, val+mcMAXt)]+1])
                labels.append('')
                lbempty.append(True)
                mcMAXt -=1
    
    
            if carbamido:
                for val, pep in enumerate(peplist):
                    if pep.count('C')>0:
                        c_locus = [i for i, letter in enumerate(pep) if letter == 'C']
                        labels[val] += f'{len(c_locus)}xCarbamidomethyl [{"; ".join(["C"+str(x+1) for x in c_locus])}]'
                        lbempty[val] = False
    
    
            if tmtlab:
                for val,pep in enumerate(peplist):
                    if len([i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))])>0:
                        tmt = [i for i, letter in enumerate(pep) if ((letter == 'K') and(i<len(pep)-1))]
                        if not lbempty[val]:
                            labels[val] += '; '
                        labels[val] += f'{len(tmt)}xTMTpro [{"; ".join(["K"+str(x+1) for x in tmt])}]'
                        lbempty[val] = False
                    if lbempty[val]:
                        labels[val] += '1xTMTpro [N-Term]'
                    else:
                        labels[val] += '; 1xTMTpro [N-Term]'
    
            pos_seq = [f"[{base+2}-{base+1+len(pep)}]" for pep in peplist]
            peplines = [[labels[val], pos_seq[val], pep, 'LysC', mcMAX-val] for val, pep in enumerate(peplist)]
    
    
            addition_dex = ['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites']
    
            pep_net = [pd.Series(pepline, index=addition_dex).transpose() for pepline in peplines]
    
            peps = peps.append(pep_net).reset_index(drop=True)


    return peps[peps['Sequence']!='']


def seq2pep(df, timesplit):
    df = df.reset_index(drop=True)
    modlist = [mod_interp(x) for x in df['Modifications']]

    timesplit = timed_print('Modifications Parsed.', timesplit)

    mod_seqlist = [translate_mods(seq, modlist[mod_id]) for mod_id, seq in enumerate(df['Sequence'])]

    timesplit = timed_print('Mods Converted to Structure.', timesplit)

    chemlist = [seq2mol(x) for x in mod_seqlist]

    timesplit = timed_print(f'{len(modlist)} Peptides Designed.', timesplit)

    df['SMILES string'] = pd.Series(chemlist)
    return df, timesplit


def mod_interp(modstring):
    modstring = str(modstring)
    mods = {'TMTn':0, 'TMT':0, 'C':0, 'O':0, 'ACn':0, 'AC':0}
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
    return mods


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


def calc_descriptors(df, timesplit):
    targets = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']
    timesplit = timed_print(f'{df.shape[0]} Peptides Written to File.', timesplit)

    chemlist = [Chem.MolFromSmiles(x) for x in df['SMILES string']]

    timesplit = timed_print(f'{df.shape[0]} Peptides Rendered.', timesplit)

    tot = len(chemlist)
    descrs = [calc_target_descriptors(mol, val, tot) for val, mol in enumerate(chemlist)]
    dfd = pd.DataFrame(descrs)

    print()
    timesplit = timed_print(f'{dfd.size} Molecular Descriptors calculated.', timesplit)

    return pd.concat([df, dfd], axis=1).dropna(subset=targets), timesplit


def calc_target_descriptors(mol, val, tot):
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


def clf_peptides(peps, clf):
    targets = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']

    peps = peps.dropna(subset=targets)

    xs = peps[targets]
    predictions = clf.predict(xs)
    pred_prob = clf.predict_proba(xs)
    peps['prediction'] = predictions
    peps['predscore'] = pd.DataFrame(pred_prob)[1][:]
#    peps['prediction lprob'] = pred_l_prob
    return peps


def score_peptides(peps, prot, enzym, mcMAX):
    prot_residues = [aa for aa in prot]
    loci = [val+1 for val, aa in enumerate(prot)]

    resplex = pd.DataFrame({'Amino Acid':prot_residues, 'Locus':loci})
    print(resplex.shape)

    average_cleaverate={0:0.6623666,
                        1:0.284661434,
                        2:0.0499406,
                        3:0.00295705,
                        4:0.0000588416,
                        5:0.0000154542,
                        6:0.0}

    res_score = [{} for _ in range(mcMAX)]
    tot_res_score = {}
    res_peps = {}
    optim_mc = {}
    temp_score = {}
    mc_denom = sum([average_cleaverate[x] for x in range(mcMAX)])
    print(mc_denom)

    for enz in enzym.keys():
        if enzym[enz]:
            tot_res_score[enz] = [0 for aa in prot]
            res_peps[enz] = ['' for aa in prot]
            optim_mc[enz] = [0 for aa in prot]
            temp_score[enz] = ['' for aa in prot]
            for i in range(mcMAX):
                res_score[i][enz] = [0 for aa in prot]

    for pep in peps.index:
        entry = peps.loc[pep][:]
        pos = str(entry['Positions in Master Proteins'])[1:-1].split('-')
        enz_a = entry['Enzyme']
        mc_a = entry['CleaveSites']
        pepdesc = f"({str(entry['Sequence'])}{str(entry['Positions in Master Proteins'])}mc{mc_a-1}, {str(entry['Modifications'])});"
        
        for x in range(int(pos[0])-1, int(pos[1])):
            res_peps[enz_a][x] += pepdesc
            temp_score[enz_a][x] += f"({str(mc_a)}, {str(entry['predscore'])});"

    for enz in enzym.keys():
        if enzym[enz]:
            score_breakdown_raw = [rep.split(';')[:-1] for rep in temp_score[enz_a]]
            score_breakdown = [[[float(b) for b in g[1:-1].split(', ')] for g in i] for i in score_breakdown_raw]
#            mc_locus = [0 for aa in prot]
            for i in range(mcMAX):
                score_analysis = [pd.DataFrame(c)[pd.DataFrame(c)[0]==i+1][1].mean() for c in score_breakdown]
                res_score[i][enz] = score_analysis

    for enz in enzym.keys():
        if enzym[enz]:
            resplex[f"{enz}: Peptides"] = res_peps[enz_a]
            ideal_mc = [0 for aa in prot]
            topscore_mc = [0 for aa in prot]
            global_score = [0 for aa in prot]
            for mc in range(mcMAX):
                resplex[f"{enz}: Score, mc{mc}"] = [x*100 for x in res_score[mc][enz]]
                ideal_mc = [mc if (x>=topscore_mc[val]) else ideal_mc[val] for val, x in enumerate(res_score[mc][enz])]
                topscore_mc = [x if (x>=topscore_mc[val]) else topscore_mc[val] for val, x in enumerate(res_score[mc][enz])]
                global_score = [global_score[val]+(x*(average_cleaverate[mc]/mc_denom)) for val,x in enumerate(res_score[mc][enz])]
            resplex[f"{enz}: Ideal MC"] = ideal_mc
            resplex[f"{enz}: Global Score"] = [x*100 for x in global_score]

    return resplex









timesplit = timer()

MODELpath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\'
MODEL = 'bESI_v2_3_1.pkl'

#ACCESSION = 'P10275'
ACCESSION = 'Q8TE02'
TMT = False
ALKYL = 'IAA'
enzym = {'Trypsin':False, 'GluC':False, 'ArgC':False, 'ProA':False, 'LysC':False, 'Chymotrypsin':False, 'Pep':True}
clMAX = 2

print(f'Residue Analysis of: {ACCESSION}\n')

clf = read_model(MODELpath+MODEL)
timesplit = timed_print('Loaded Learned Model.', timesplit)


protein = uniprot_call(ACCESSION)
timesplit = timed_print('Accessed Protein Information from UniProt.', timesplit)


peptides = silico_digest(protein['fasta'], TMT, ALKYL, enzym, clMAX)



timesplit = timed_print('Modelled Protein Digest.', timesplit)

pepsmiles, timesplit = seq2pep(peptides, timesplit)


physpeps, timesplit = calc_descriptors(pepsmiles, timesplit)



sortedpeps = clf_peptides(physpeps, clf)
timesplit = timed_print('New Peptides Classified.', timesplit)

sortedpeps.to_csv('pulls\\ELPpeps.csv', index=False)


residue_scores = score_peptides(sortedpeps, protein['fasta'], enzym, clMAX)
residue_scores.to_csv('pulls\\ELPres.csv', index=False)




