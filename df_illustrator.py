# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 12:58:53 2023

@author: HUFFMP
"""

import pandas as pd
from time import sleep as slp
from time import time
import sys
import requests
import re
from requests.adapters import HTTPAdapter, Retry

start = time()

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

print('Which columns would you like added?')
yn_gene = input('Gene Name (y/n): ').lower()=='y'
yn_allgene = input('All Gene Names (y/n): ').lower()=='y'
yn_prot = input('Full Protein Name (y/n): ').lower()=='y'
yn_allprot = input('All Protein Names (y/n): ').lower()=='y'
yn_len = input('Protein Length (y/n): ').lower()=='y'


lst_gene = []
lst_allgene = []
lst_prot = []
lst_allprot = []
lst_len = []


count = 0
print('Note: Sometimes this can take awhile.')
print(f"Estimated total time: {round((time()-start)+0.214*len(df.index), 2)} seconds")
for acc in df[df.columns[int(ab_select)-1]]:
    EXISTS = True
    acc = acc.split('-')[0].split(';')[0].split(' ')[0]
    url = f'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+accession:{str(acc)}&format=tsv'
#    print(url)
    for batch, total in get_batch(url):
        dict1 = {}
        for key, val in enumerate(batch.text.splitlines()[0].split('\t')):
            if len(batch.text.splitlines()) > 1:
                dict1[val] = batch.text.splitlines()[1].split('\t')[key]
            else:
                EXISTS = False
#        print(dict1)
    if EXISTS:
        if yn_gene:
            lst_gene.append(dict1['Gene Names'].split(' ')[0])
        if yn_allgene:
            lst_allgene.append('; '.join(dict1['Gene Names'].split(' ')))
        if yn_prot:
            lst_prot.append(dict1['Protein names'].split('(')[0])
        if yn_allprot:
            lst_allprot.append(dict1['Protein names'])
        if yn_len:
            lst_len.append(dict1['Length'])
    else:
        if yn_gene:
            lst_gene.append('')
        if yn_allgene:
            lst_allgene.append('')
        if yn_prot:
            lst_prot.append('')
        if yn_allprot:
            lst_allprot.append('')
        if yn_len:
            lst_len.append('')
    print('.',end='')

print()
if yn_gene:
    df['Gene Name'] = lst_gene
if yn_allgene:
    df['All Gene Names'] = lst_allgene
if yn_prot:
    df['Full Protein Name'] = lst_prot
if yn_allprot:
    df['All Protein Names'] = lst_allprot
if yn_len:
    df['Protein Length'] = lst_len


dataname = '.'.join(datasource.split('.')[:-1])
df.to_excel(f"{dataname}_illustrated.xlsx", index=False)
print(f"{time()-start} seconds.")



