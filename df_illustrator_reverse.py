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

ab_select = int(input('Enter the index of the GENE NAME column for analysis.\n'+
                  'Example: if column 4 is desired, input should be: 4\nIndex (-1 to exit): '))-1
if ab_select==-1:
    stop()

print('Which columns would you like added?')
yn_gene = input('SWISSPROT accession (y/n): ').lower()=='y'


lst_gene = []


count = 0
print(f"Estimated total time: {(time()-start)+0.185*len(df.index)} seconds")
print(df[df.columns[int(ab_select)]])
for acc in df[df.columns[int(ab_select)]]:
    acc = acc.split('-')[0].split(';')[0].split(' ')[0]
    url = f'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+gene:{str(acc)}&format=tsv'
#    print(url)
    for batch, total in get_batch(url):
        valid = True
        dict1 = {}
        if len(batch.text.splitlines())>1:
            for key, val in enumerate(batch.text.splitlines()[0].split('\t')):
                dict1[val] = batch.text.splitlines()[1].split('\t')[key]
        else:
            valid = False
#        print(dict1)
    if yn_gene and valid:
        lst_gene.append(dict1['Entry'].split(' ')[0])
    elif yn_gene:
        lst_gene.append('')
    count += 1
#    print(f"{count}/{len(df.index)}")

if yn_gene:
    df['accession'] = lst_gene


datasite = '\\'.join(datasource.split('\\')[:-1])
datacode = datasource.split('\\')[-1].split('.')[0]
df.to_csv(f"{datasite}\\{datacode}_illustrated.csv", index=False)


print(f"{time()-start} seconds.")



