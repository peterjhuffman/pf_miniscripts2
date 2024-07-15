# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:04:04 2022

@author: HUFFMP
"""
# \\\\amer.pfizer.com\\pfizerfiles\\Research\LAJ\\oru-tcb\\ChemBio\\MS_Analyses_PD\\SKP2\\20221205_E3OverExpression_global_PH\\volcanos\\btrc.csv
import pandas as pd
from time import sleep as slp
import sys
import requests
import re
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

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
    print('Exiting spotifyre.py. Thanks for dropping by.')
    slp(2)
    sys.exit()

uniprotID = ''
datasource = ''

datasource = input('DATA LOCATION (x to exit): ')
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

# UNIPROT CALL -----------------------------------------------------------------------------
#uniprotID = input('UNIPROT ID (x to exit): ')
#if uniprotID.lower()=='x':
#    stop()
#
#
#url = f'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query={uniprotID}%20AND%20%28reviewed%3Atrue%29&size=500'
#interactions = {}
#for batch, total in get_batch(url):
#    for line in batch.text.splitlines()[1:]:
#        primaryAccession, interactsWith = line.split('\t')
#        interactions[primaryAccession] = len(interactsWith.split(';')) if interactsWith else 0
#    print(f'{len(interactions)} / {total}', end = ' interactions.    ')
#
#print('Interactome Loaded.')
#itxome = list(interactions.keys())
# UNIPROT CALL -----------------------------------------------------------------------------

subsource = input('KNOWN SUBSTRATES (x to exit): ')
if datasource.lower()=='x':
    stop()
#itxome = list(pd.read_csv(subsource, sep='\t')['Gene Symbol(Substrate)'])
itxome = list(pd.read_csv(subsource, sep='\t')['SUBGENE'])
print(itxome)

subsource = input('PREDICTED SUBSTRATES (x to exit): ')
if datasource.lower()=='x':
    stop()
#itxome = list(pd.read_csv(subsource, sep='\t')['Gene Symbol(Substrate)'])
itxome_p = list(pd.read_csv(subsource, sep='\t')['Gene Symbol(Substrate)'])
print(itxome_p)


df['color'] = ''
df['label'] = ''
df['accession2'] = ''
df['accession3'] = ''
df['abundance'] = ''

#df['accession2'] = [x[0] for x in df['Accession'].str.split('-')]
df['accession3'] = [x[1] for x in df['Protein ID'].str.split('|')]
df['accession2'] = [x[0] for x in df['accession3'].str.split('-')]

THRESH=.10
CUTOFF = pd.Series(df[df.columns[-23:-5]].to_numpy().flatten()).quantile(THRESH)
df['abundance'] = (df[df.columns[-23:-5]].mean(axis=1)>CUTOFF).astype(int)

mask_p = df['the -LOG(P-value)'] > 1.301
diff_std = df['Difference'].std()*3

plus_mask = df['Difference'] < (0-diff_std)
minus_mask = df['Difference'] > (diff_std)
either_mask = (df['Difference'] > (diff_std))|(df['Difference'] < (0-diff_std))

abundance_mask = df['abundance']>0

#mask_p = df['Significant'] == '+'
#plus_mask = df['Difference'] < 0
#minus_mask = df['Difference'] > 0

df.loc[mask_p&plus_mask&abundance_mask, 'label'] = df[mask_p&plus_mask]['Gene Symbol']
df.loc[mask_p&plus_mask, 'color'] = 'shift'

df.loc[mask_p&minus_mask&abundance_mask, 'label'] = df[mask_p&minus_mask]['Gene Symbol']
df.loc[mask_p&minus_mask, 'color'] = 'shift'


#df.loc[df['Gene'].isin(itxome), 'label'] = df[df['accession2'].isin(itxome)]['Gene']

df.loc[df['Gene Symbol'].isin(itxome_p)&mask_p&either_mask, 'color'] = 'itx_p'
print(f"{df.loc[df['Gene Symbol'].isin(itxome_p)&mask_p&either_mask].shape[0]} significant predicted substrates identified through MS.")

df.loc[df['Gene Symbol'].isin(itxome)&mask_p&either_mask, 'color'] = 'itx'
print(f"{df.loc[df['Gene Symbol'].isin(itxome)&mask_p&either_mask].shape[0]} significant substrates identified through MS.")




newfile = f"{'.'.join(datasource.split('.')[:-1])}_itxome2.csv"

print('Saving file.')
df.to_csv(newfile, index=False)