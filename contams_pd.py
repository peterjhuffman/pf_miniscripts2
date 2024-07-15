# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 12:10:47 2023

@author: HUFFMP
"""


import pandas as pd
import pathlib

# json_loc = sys.argv[1]

# with open(json_loc, encoding="utf-8") as g:
#     lines = json.load(g)

# infoloc1 = lines['Tables'][0]['DataFile']
# infoloc2 = lines['Tables'][1]['DataFile']
# jobname = lines['ResultFilePath'].split('.')[0].split('/')[-1]

# df = pd.read_csv(infoloc1, sep='\t')
# sdf = pd.read_csv(infoloc2, sep='\t')


filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\contaminants2.fasta'
f = open(filepath, "r")
lines = f.readlines()
f.close()

toplines = []
for line in lines:
    if line.count('>')>0:
        toplines.append(line)

ids = []
for line in toplines:
    ids.append(line[1:7])


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
    input('Restart and Try again.')
    print(nothing)
print('File read successfully.')
v_target = filename

df['accession3'] = 0
df.loc[df['Accession'].isin(ids), 'accession3'] = 1


df.loc[df['Description'].str.count('ENSEMBL')>0, 'accession3'] = 1
df.loc[df['Description'].str.count('contam')>0, 'accession3'] = 1
df.loc[df['Description'].str.count('taurus')>0, 'accession3'] = 1
df.loc[df['Description'].str.count('SWISS')>0, 'accession3'] = 1

df = df[df['accession3']==0].reset_index(drop=True)
del df['accession3']
df.to_excel(v_target[:-len(fileSuffix)]+'_pure.xlsx')
input('Done. press any key to exit')