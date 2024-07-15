# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:08:37 2024

@author: HUFFMP
"""

import numpy as np
import pandas as pd

anno = pd.read_csv('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan slices/model anno.csv')

nsclc = [x for x in anno[anno['cancer_type']=='Non-Small Cell Lung Carcinoma']['model_name'].values]

tolu = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu taa/Complete TAA List for NSCLC_ill.xlsx')
# toltargets = [x for x in tolu['accession'].values]

# proc = pd.read_csv('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan slices/source/Protein_matrix_averaged_20221214.tsv', sep='\t', low_memory=False)


# del proc['Unnamed: 1']
# proc=proc.T
# proc.index = proc[0]
# del proc[0]
# del proc[1]


# proc_c = proc[proc['uniprot_id'].isin(nsclc)].T

# proc_t = proc_c[proc_c.index.isin(toltargets+['uniprot_id'])]







# proc_t = proc[proc.index.isin([x for x in tolu['Gene Symbol'].values])]

# rna = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan slices/source/procan rnaseq tpm.xlsx')

# rna_t = rna[rna['symbol'].isin([x for x in tolu['Gene Symbol'].values])]
rna_c = rna_t[['symbol']+list(set(nsclc).intersection(set([x for x in rna_t.columns])))]





























