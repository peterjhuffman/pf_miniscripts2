# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:20:28 2023

@author: HUFFMP
"""

# Imports
from sklearn.datasets import make_blobs
from sklearn.model_selection import train_test_split
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import confusion_matrix
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.decomposition import PCA
from time import time
from random import sample

frac_test_split = 0.33
random_seed = 14

def scale(df, keycards):
    x = df[keycards].values.reshape(-1, 1)
    y = df['Acquired'].values

    bal = 1/(y.sum()/len(y))
#    print(bal)

    # Generate data
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=frac_test_split, random_state=random_seed)
    
#    vx, vy = 0, 1
    
#    # Generate scatter plot for training data 
#    plt.scatter(X_train[:,vx], X_train[:,vy])
#    plt.title('Linearly separable data')
#    plt.xlabel(keycards[vx])
#    plt.ylabel(keycards[vy])
#    plt.show()
#    #plt.clf()
#    
#    plt.scatter(X_train[y_train==0][:,vx], X_train[y_train==0][:,vy], color='green')
#    plt.scatter(X_train[y_train==1][:,vx], X_train[y_train==1][:,vy], color='purple')
#    plt.title('Linearly separable data')
#    plt.xlabel(keycards[vx])
#    plt.ylabel(keycards[vy])
#    plt.show()
    
    
    
    # Initialize SVM classifier
#    print(keycards)
    wetweight=1
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto')
    clf = clf.fit(X_train, y_train)
    
    
    # Predict the test set
    predictions = clf.predict(X_test)
    
    print('',end='-')
    # Generate confusion matrix
    matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)
    
#    print(matrix)
    acc_tot = round((matrix[1][1]+matrix[0][0])/matrix.sum()*100, 1)
    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)

#    print(f'Total Accuracy: {acc_tot}%')
#    print(f'Target Accuracy: {acc_targ}%')

    return acc_tot, acc_targ


def scale_add(df, descrs, newkeys, max_iter):
    x = df[newkeys].values
    y = df['Acquired'].values

    bal = 1/(y.sum()/len(y))
    #    print(bal)

    # Generate data
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=frac_test_split, random_state=random_seed)


    # Initialize SVM classifier
    wetweight=1
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto')
    clf = clf.fit(X_train, y_train)

    # Predict the test set
    predictions = clf.predict(X_test)

    # Generate confusion matrix
    matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)

    acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)

#    print(newkeys)
#    print(f'Total Accuracy: {acc_tot}%')
#    print(f'Target Accuracy: {acc_targ}%')
#    print()

#    for key in descrs.sort_values('acc rank')['keys'].values:
#    for key in descrs.sample(frac=1)['keys'].values:
    for key in descrs:
        black_keys = newkeys+[key]
        x = df[black_keys].values
        y = df['Acquired'].values
    
        bal = 1/(y.sum()/len(y))
        #    print(bal)
    
        # Generate data
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=frac_test_split, random_state=random_seed)
    
    
        # Initialize SVM classifier
        wetweight=1.0
        weights = {0:1.0, 1:bal*wetweight}
        clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto')
        clf = clf.fit(X_train, y_train)
    
        # Predict the test set
        predictions = clf.predict(X_test)
    
        # Generate confusion matrix
        matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)
    
        bacc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
        bacc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)

        if ((bacc_tot + bacc_targ) > (acc_tot + acc_targ)):
            newkeys.append(key)
            max_iter -= 1
            if max_iter > 0 :
                newkeys = scale_add(df, descrs, newkeys, max_iter)
            else:
                return newkeys
            break
    return newkeys


def scale_maxer(df, descrs, newkeys, max_iter, repeat):
    acc_rank = 0
    g = repeat


    while repeat >= 1:
        descx = descrs.sample(frac=1)['keys'].values
        newkeys = scale_add(df, descx, newkeys, max_iter)
    
        newkeys = list(set(newkeys))
        print(newkeys)
        x = df[newkeys].values
        y = df['Acquired'].values
        
        bal = 1/(y.sum()/len(y))
        #    print(bal)
        
        # Generate data
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=frac_test_split, random_state=random_seed)
        
        
        # Initialize SVM classifier
        wetweight=1.0
        weights = {0:1.0, 1:bal*wetweight}
        clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto')
        clf = clf.fit(X_train, y_train)
        
        
        # Predict the test set
        predictions = clf.predict(X_test)
        
        
        # Generate confusion matrix
        matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)
        acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
        acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)
        print(f'Total Accuracy: {acc_tot}%')
        print(f'Target Accuracy: {acc_targ}%')
        print()
        
        if acc_tot + acc_targ > acc_rank:
            keyMAX = newkeys
            loca = repeat
            acc_rank = acc_tot + acc_targ

        repeat -=1


    print(f'Target selected on round {g-loca+1}')
    return keyMAX


def read_file(datasource):
    fdtype = datasource[-3:]
    if fdtype == 'csv':
        df = pd.read_csv(datasource)
    elif fdtype == 'lsx':
        df = pd.read_excel(datasource)
    return df.dropna()


def df_compile(list_files):
    print(list_files[0][104:-10], end='} ')
    df = read_file(list_files[0])
    for file in list_files[1:]:
        print(file[104:-10], end='} ')
        df_a = read_file(file)
        df = df.append(df_a).reset_index(drop=True)
    return df

def pca_tab(df):
    m = Chem.MolFromSmiles('C')
    vals = Descriptors.CalcMolDescriptors(m)
    keycards = list(vals.keys())
    x = df[keycards].astype(np.float64)
    
    pca = PCA(n_components=3)
    
    pCones = pca.fit(x.values)

    principalDf = pd.DataFrame(pCones, columns = ['PC1', 'PC2', 'PC3'])
    df['PC1'] = principalDf['PC1']/max(principalDf['PC1'])
    df['PC2'] = principalDf['PC2']/max(principalDf['PC2'])
    df['PC3'] = principalDf['PC3']/max(principalDf['PC3'])

    return df

def SVM_test(df, newkeys):
    print()
    newkeys = list(set(newkeys))
    print(newkeys)
    
    print(str(len(newkeys))+' descriptors analyzed')
    x = df[newkeys].values
    y = df['Acquired'].values
    
    bal = 1/(y.sum()/len(y))
    #    print(bal)
    
    # Generate data
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=frac_test_split, random_state=random_seed)
    
    
    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto')
    clf = clf.fit(X_train, y_train)
    
    
    # Predict the test set
    predictions = clf.predict(X_test)
    
    
    # Generate confusion matrix
    matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)
    
    print(matrix)
    acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)
    
    print(f'Total Accuracy: {acc_tot}%')
    print(f'Target Accuracy: {acc_targ}%')
    print(f'Accuracy Rank: {acc_targ+acc_tot}')

    return acc_tot, acc_targ


def random_key(keys, sample_v):
    return sample([x for x in keys], sample_v)


RUNS = 12
MAXDF_SIZE = 30
MINDF_SIZE = 5
VOLUME = 10
REPS = 1
REAL_REPS = 2

clock = time()

#randpic_check = sample([x for x in range(1, 178)], MAXDF_SIZE)
randpic_check = [42, 43, 44, 123, 124, 125, 82, 127, 83, 128, 164, 126, 41, 84, 129, 130, 45]

path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\pepprtmt\\desc\\'
files = [path + f'PFEC_PEPPRpep_{x}_desc.xlsx' for x in randpic_check]
dfmax = df_compile(files)
print(dfmax.shape)

m = Chem.MolFromSmiles('C')
img = Chem.Draw.MolToImage(m)
vals = Descriptors.CalcMolDescriptors(m)
keycards = list(vals.keys())
descrs = pd.DataFrame()
descrs['keys'] = keycards

#key_acc = [scale(dfmax, keycard) for keycard in keycards]
#descrs['total acc'] = [key[0] for key in key_acc]
#descrs['target acc'] = [key[1] for key in key_acc]
#descrs['acc rank'] = descrs['target acc'] + descrs['total acc']
#key_sample = descrs.sort_values('acc rank')['keys'].values[-60:]


acc_rank = 0
for x in range(RUNS):
    print()
    print()
    print()
    randpic = sample(randpic_check, 2)
    files = [path + f'PFEC_PEPPRpep_{x}_desc.xlsx' for x in randpic]
    df = df_compile(files)
    print(df.shape)
    key_acc = [scale(df, keycard) for keycard in keycards]
    descrs['total acc'] = [key[0] for key in key_acc]
    descrs['target acc'] = [key[1] for key in key_acc]
    descrs['acc rank'] = descrs['target acc'] + descrs['total acc']
#    randpic = sample([x for x in range(1, 178)], MINDF_SIZE)

    #descrs.to_csv('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\key_rank.csv',index=False)
    #descrs = pd.read_csv('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\key_rank.csv')
    
    #newkeys = descrs3['keys'].values
#    newkeys = ['AvgIpc', 'SlogP_VSA10', 'EState_VSA8']
    #newkeys = ['Chi3v', 'Chi1v', 'FpDensityMorgan3', 'VSA_EState10']
    #newkeys = ['EState_VSA4', 'PEOE_VSA4', 'fr_bicyclic', 'SlogP_VSA8', 'fr_quatN', 'MinEStateIndex', 'FpDensityMorgan2']
    #newkeys = ['EState_VSA11', 'fr_unbrch_alkane', 'SlogP_VSA9', 'Chi4v', 'fr_para_hydroxylation', 'fr_priamide', 'fr_phenol_noOrthoHbond', 'SlogP_VSA8', 'Chi4n', 'fr_Ar_OH', 'MinEStateIndex', 'fr_Ar_N', 'fr_guanido']

    for descrval in [x for x in descrs.sort_values('acc rank')['keys'].values[-REAL_REPS:]]:
        startkey = [descrval]
    
        print(f'Startkey: {startkey[0]}')
        newkeys = scale_maxer(df, descrs, startkey, VOLUME, REPS)
    #    newkeys = scale_maxer(df, descrs, newkeys, 8, 3)
        print('===============================================================')
        
        acc_tot, acc_targ = SVM_test(dfmax, newkeys)
        if acc_rank < acc_tot+acc_targ:
            glob_keys = newkeys
            acc_rank = acc_tot+acc_targ
        print()




print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(glob_keys)
SVM_test(dfmax, glob_keys)
print()
print(f'Total runtime: {round(time()-clock, 1)} sec')
