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


# Setting up URL grab settings for UNIPROT call
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))



def timer_quant(sec):
    """
    Splits a value in seconds into hours, minutes and seconds.
    Receives:
        sec : total seconds
    Returns:
        n_hr : number of hours
        n_min : number of minutes
        n_sec : number of seconds
    """
    nclock = round(sec, 1)
    n_hr = int(nclock/3600)
    n_min = int((nclock-(3600*n_hr))/60)
    n_sec = nclock - (3600*n_hr) - 60*n_min
    return (n_hr, n_min, n_sec)


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


def get_next_link(headers):
    """
    Architecture for calling uniprot. From https://www.uniprot.org/help/api_queries
    Receives:
        not sure
    Returns:
        not sure
    """
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    """
    Architecture for calling uniprot. From https://www.uniprot.org/help/api_queries
    Receives:
        not sure
    Returns:
        not sure
    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def stop():
    """
    Ends program.
    Receives:
        
    Returns:
        
    """
    print('Exiting protein_ESI.py. Thanks for dropping by.')
    slp(2)
    sys.exit()


def read_model(filename):
    """
    Reads the classifier model saved as a pickle file.
    Receives:
        filename : the full file location of the classifier model file.
    Returns:
        clf_R : the classifier
    """
    with open(filename, 'rb') as fid:
        clf_R = pickle.load(fid)
    return clf_R


def uniprot_call(acc):
    """
    Calls uniprot for the accession file.
    Receives:
        acc : a swissprot accession number like : "P10275" for AR
    Returns:
        d1 : a dictionary with: 
            ['header'] : uniprot header line
            ['fasta'] : protein fasta
    """
    url = f'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+(organism_id:9606)+AND+accession:{str(acc)}&format=fasta'
    for batch, total in get_batch(url):
        d1 = {}
        d1['header'] = batch.text.splitlines()[0]
        d1['fasta'] = ''.join(batch.text.splitlines()[1:])
    return d1


def silico_digest(stat_prot, tmtlab, alkyl, enzym, mcMAX):
    """
    This is the in-silico digestion algorithm. turns a fasta file into a DataFrame of peptides depending
    on a number of received characteristics. Currently enabled for TMT, carbimidomethyl, alkylation.
    Hopefully will add more in the future.
    Receives:
        stat_prot : the fasta file of the protein of interest. text
        tmtlab : a boolean - whether the new peptides should be TMT labeled or not
        alkyl : the reagent used to alkylate the peptides. Currently only 'IAA' and 'CAA' are supported.
        enzym : a dictionary of possible enzymes for digestion. current enzymes supported:
            trypsin, chymotrpysin, gluc, argc, proa. lysC is not currently functional.
        mcMAX : number of cleave site residues allowed. Minimum is 1 (meaning a peptide with one cleavesite, no missed cleavages).
            I recommend using only 1 or 2 right now - the model was trained on peptides with either 0 or 1 MCs.
            2+ MCs get pretty big and cause overflow problems, plus are <5% of all peptides anyway
    Returns:
        peps: A dataframe of information on each peptide. Dataframe contains info on peptide:
            Modifications : in standard PD output
            Positions in Master Proteins : also in standard PD output
            Sequence : raw text sequence
            Enzyme : enzyme used
            CleaveSites : number of cleavesites in peptides (MC + 1)
    """
    carbamido = ((alkyl=='IAA') or (alkyl=='CAA'))
    mcMAXt = mcMAX

    peps = pd.DataFrame(columns=['Modifications', 'Positions in Master Proteins', 'Sequence', 'Enzyme', 'CleaveSites'])

    prot=stat_prot

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
            peps = peps[peps['Sequence']!='']


        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)

    prot=stat_prot
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
            peps = peps[peps['Sequence']!='']

        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)


    prot=stat_prot
    if enzym['ArgC']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = prot[j+1:].find('R')
            if i==-1:
                i = prot[j+1:].find('R')
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
            peps = peps[peps['Sequence']!='']

        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)


    prot=stat_prot
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
            peps = peps[peps['Sequence']!='']

        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)


    prot=stat_prot
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
            peps = peps[peps['Sequence']!='']

        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)


    prot=stat_prot
    if enzym['LysC']:
        bases = [-1]
        i, j = 0, -1
        while i !=-1:
            i = prot[j+1:].find('K')
            if i==-1:
                i = prot[j+1:].find('K')
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
            peps = peps[peps['Sequence']!='']

        peps=pd.concat([peps.iloc[:len(peps.index)-mcMAX],pd.DataFrame(peps.iloc[len(peps.index)-1]).transpose()]).reset_index(drop=True)



    return peps[peps['Sequence']!='']


def seq2pep(df, timesplit):
    """
    Receives a DataFrame of peptides and adds a column of peptide in SMILES format.
    Receives:
        df : dataframe of peptides, outputted from a silico digest
        timesplit : time of last message
    Returns:
        df : dataframe of peptides, column containing SMILES added
        timesplit : time of last message
    """
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
    """
    Translates a PD output modification string into a dictionary containing information on modifications.
    Receives:
        modstring : string output from PD 'Modifications'
    Returns:
        mods: dict object containing info on TMT, Acetyl, alkyl and oxidation
    """
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
    """
    Receives an amino acid sequence and returns a SMILES string.
    Receives:
        seq : string of amino acids
    Returns:
        smile : string of SMILES code for peptide
    """
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
    """
    Modifies string of amino acids to contain information on modifications.
    Receives:
        seq : sequence of canonical amino acids
        mods: dict object containing info on TMT, Acetyl, alkyl and oxidation
    Returns:
        seq : modified sequence with non-canonical modified aminos substituted in.
            c : carbamidomethyl
            o : oxidize methionine
            m : n-terminal TMT
            y : n-terminal acetyl
            t : lysine tmt
            a : acetyllysine
    """
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


def mol_define(seq):
    if len(seq)<900:
        return Chem.MolFromSmiles(seq)
    else:
        return 'x'


def calc_descriptors(df, timesplit):
    """
    Calculates 12 molecular descriptors for a dataframe of peptides based on SMILES sequences.
    These molecular descriptors were selected with a variable selection genetic-ish algorithm.
    Receives:
        df : dataframe of peptides containing a column 'SMILES string' of peptide SMILES representation
        timesplit: time since last message
    Returns:
        df : dataframe of peptides containing 12 molecular descriptors. 
        timesplit : time since last message
    """
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


def calc_target_descriptors(mol, val, tot):
    """
    Calculates 12 molecular descriptors using the RDKit library.
    Receives:
        mol : molecule in rdkit.Chem form
        val : number of peptide in database for message call
        tot : total number of peptides in database for message call
    Returns:
        dx : dictionary of molecular descriptors of mol
    """
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


def clf_peptides(peps, clf):
    """
    Classifies peptides according to learned model and molecular descriptors.
    Receives:
        peps : dataframe of peptides, containing descriptor values.
        clf : loaded classifier
    Returns:
        peps : dataframe of peptides, containing added columns:
            ['prediction'] : 0 or 1, prediction made by classifier
            ['predscore'] : prediction probability made by classifier
    """
    targets = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc',
               'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']

    xs = peps[targets]
    predictions = clf.predict(xs)
#    pred_prob = clf.predict_proba(xs)
    peps['prediction'] = predictions
#    peps['predscore'] = pd.DataFrame(pred_prob)[1][:] # ------------------------------------------------------------------------------------------------------------
#    peps['prediction lprob'] = pred_l_prob

    for ind in list(peps.loc[(peps['fr_Nhpyrrole']==0)&
                             (peps['fr_guanido']==0)&
                             (peps['BCUT2D_CHGLO']==0)&
                             (peps['MinEStateIndex']==0)&
                             (peps['Chi4v']==0)&
                             (peps['fr_ether']==0)&
                             (peps['AvgIpc']==0)&
                             (peps['fr_unbrch_alkane']==0)&
                             (peps['BCUT2D_MWLOW']==0)&
                             (peps['BCUT2D_MRHI']==0)&
                             (peps['BCUT2D_LOGPHI']==0)&
                             (peps['FpDensityMorgan1']==0)].index):
        peps.loc[ind]['prediction'] = 0
        peps.loc[ind]['predscore'] = 0


    return peps


def score_peptides(peps, prot, enzym, mcMAX):
    """
    Scores individual amino acids based on peptide scoring.
    Receives:
        peps : dataframe of peptides, containing peptide scoring
        prot : a dictionary with: 
            ['header'] : uniprot header line
            ['fasta'] : protein fasta
        enzym : a dictionary of possible enzymes for digestion. current enzymes supported:
            trypsin, chymotrpysin, gluc, argc, proa. lysC is not currently functional.
        mcMAX : number of cleave site residues allowed. Minimum is 1 (meaning a peptide with one cleavesite, no missed cleavages).
    Returns:
        resplex : a dataframe of amino acids with scores for each digest. + scores for each MC within each digest.
    """

    pred_id = 'prediction'

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
            temp_score[enz_a][x] += f"({str(mc_a)}, {str(entry[pred_id])});"


    for enz in enzym.keys():
        if enzym[enz]:
            print(enz)
            score_breakdown_raw = [rep.split(';')[:-1] for rep in temp_score[enz]]
            score_breakdown = [[[float(b) for b in g[1:-1].split(', ')] for g in i] for i in score_breakdown_raw]
#            mc_locus = [0 for aa in prot]
            for i in range(mcMAX):
                score_analysis = [pd.DataFrame(c)[pd.DataFrame(c)[0]==i+1][1].mean() for c in score_breakdown]
                res_score[i][enz] = score_analysis

    for enz in enzym.keys():
        if enzym[enz]:
            resplex[f"{enz}: Peptides"] = res_peps[enz]
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




#base = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI tmtrtsRG-test\\AR5degrader_22RV1_RTS-peps.xlsx')
#base = base[base['Master Protein Accessions'].str.count(';')==0]
#ACCs = base['Master Protein Accessions'].unique()[165:]
#print(ACCs)
ACCs = ['Q15047']



T,F = True, False
timesplit = timer()

# where model is stored
MODELpath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\'
#model name
MODEL = 'bESI_v2_5_0p.pkl'
# protein of interest
#ACCESSION = 'P10275'
#ACCESSION = 'Q8TE02'

# info on modifications
TMT = False
ALKYL = 'IAA'

#enzymes used
enzym = {'Trypsin':T, 'GluC':F, 'ArgC':F, 'ProA':F, 'LysC':F, 'Chymotrypsin':F}
#maximum number of cleave sites
clMAX = 2 #recommended 2 or below. 3+ gives more information, but model is less accurate


posi = []
neg = []
acc_game = []

for ACCESSION in ACCs:
    print(f'Residue Analysis of: {ACCESSION}')
    print(f'TMT labeled: {TMT}')
    print(f'Enzymes: {enzym}')
    print(f'Missed Cleavage Max: {clMAX-1}\n')
    
    #loads model
    clf = read_model(MODELpath+MODEL)
    timesplit = timed_print('Loaded Learned Model.', timesplit)
    
    # obtains protein info from uniprot
    protein = uniprot_call(ACCESSION)
    timesplit = timed_print('Accessed Protein Information from UniProt.', timesplit)
    
#    digests protein using written settings
    peptides = silico_digest(protein['fasta'], TMT, ALKYL, enzym, clMAX)
#    peptides.to_csv('peptide_test.csv',index=False)
    peptides = pd.read_csv('peptide_test.csv')
    timesplit = timed_print('Modelled Protein Digest.', timesplit)
    
    # creates SMILES from peptide info
    pepsmiles, timesplit = seq2pep(peptides, timesplit)
    
    # calculates descriptors of peptides
    physpeps, timesplit = calc_descriptors(pepsmiles, timesplit)
    
    tmt_code = {True:'_TMT', False:''}
    # scores peptide according to classifier
    sortedpeps = clf_peptides(physpeps, clf)
    timesplit = timed_print('New Peptides Classified.', timesplit)
#    sortedpeps.to_csv(f'pulls\\custom_res.csv', index=False)

#    sortedpeps.to_csv(f'pulls\\{ACCESSION}{tmt_code[TMT]}peps.csv', index=False)
#    sortedpeps = pd.read_csv(f'pulls\\{ACCESSION}{tmt_code[TMT]}peps.csv')
#    scores residues according to scored peptides
    residue_scores = score_peptides(sortedpeps, protein['fasta'], enzym, clMAX)
    residue_scores.to_csv(f'pulls\\{ACCESSION}{tmt_code[TMT]}res.csv', index=False)
#
#    dfa = pd.DataFrame(columns=[ACCESSION])
#
#    for acc in [ACCESSION]:
#        protein = uniprot_call(acc)
#        sliced = base[base['Master Protein Accessions']==acc]
#        psms = [0 for aa in protein['fasta']]
#    
#        for pep in sliced.index:
#            entry = sliced.loc[pep][:]
#            pos = str(entry['Positions in Master Proteins'])[1:-1].split('[')[1].split('-')
#            for p in range(int(pos[0])-1, int(pos[1])):
#                psms[p] += int(entry['# PSMs'])
#        dfa[acc] = pd.Series(psms)
#
#    dfa['score'] = residue_scores['Trypsin: Global Score']
#
#    dfa_0 = dfa[dfa[acc]==0]
#    dfa_1 = dfa[dfa[acc]!=0]
#
#    dfa_00 = dfa_0[dfa_0['score']==0]
#    dfa_01 = dfa_0[dfa_0['score']!=0]
#    dfa_10 = dfa_1[dfa_1['score']==0]
#    dfa_11 = dfa_1[dfa_1['score']!=0]
#
#    print(dfa_11.shape[0]/(dfa_01.shape[0]+dfa_11.shape[0]))
#    print(dfa_00.shape[0]/(dfa_10.shape[0]+dfa_00.shape[0]))
#
#    pospos = dfa_11.shape[0]/(dfa_01.shape[0]+dfa_11.shape[0])
#    negneg = dfa_00.shape[0]/(dfa_10.shape[0]+dfa_00.shape[0])
#
#    posi.append(pospos)
#    neg.append(negneg)
#    acc_game.append(ACCESSION)
#
#    savegame = pd.DataFrame({'Accession':acc_game, 'Hit Odds':posi, 'Miss Odds':neg})
#    savegame.to_excel('setDB_construct.xlsx', index=False)