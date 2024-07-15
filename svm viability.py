# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:20:28 2023

@author: HUFFMP
"""

import pandas as pd 
import matplotlib.pyplot as plt
import pickle
from sklearn import svm
from sklearn.metrics import confusion_matrix


def read_model(filename):
    with open(filename, 'rb') as fid:
        clf_R = pickle.load(fid)
    return clf_R


def clf_test(df, clf, newkeys):
    x, y = df[newkeys].values, df['Acquired'].values
    predictions = clf.predict(x)
    df['predictions'] = predictions
    matrix = confusion_matrix(y, predictions, labels=clf.classes_)
    
    print(matrix)
    acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)
    
    print(f'Total Accuracy: {acc_tot}%')
    print(f'Target Accuracy: {acc_targ}%')
    print(f'Accuracy Rank: {acc_targ+acc_tot}')

    return acc_tot, acc_targ


path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI hela-test\\desc\\'

g = 60
df = pd.read_excel(path+f'PFEC_HELApep_{str(g)}_desc.xlsx').dropna()

df['len'] = [len(str(x)) for x in df['Sequence']]

df_0 = df[df['Acquired']==0]
df_1 = df[df['Acquired']==1]

plt.hist(df_0['len'], 70)
plt.show()
plt.clf()
plt.hist(df_1['len'], 20)
plt.show()
plt.clf()


i = 0
odds = []
while i < 50:
    print(i)
    df_slice = df[df['len']>i]
    df_sliced = df_slice[df_slice['len']<(i+5)]
    print(df_sliced['Acquired'].mean())
    odds.append(df_sliced['Acquired'].mean())
    i += 5


df['class'] = 0
df['class'][df['len']>10] = 1
df['class'][df['len']>20] = 0
df['class'][df['len']>35] = 1
df['class'][df['len']>45] = 0
print(df['class'].sum())

matrix = confusion_matrix(df['Acquired'], df['class'])

print(matrix)
acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)

print(f'Total Accuracy: {acc_tot}%')
print(f'Target Accuracy: {acc_targ}%')
print(f'Accuracy Rank: {acc_targ+acc_tot}')


with open(path[:-20]+'bESI_v2_0.pkl', 'rb') as fid:
    clf_R = pickle.load(fid)
clf_test(df, clf_R, ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'FpDensityMorgan1', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MRHI', 'BCUT2D_MWLOW', 'BCUT2D_LOGPHI'])
