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
import pickle

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
    dfc = []
    for file in list_files:
        print(file[104:-10], end='} ')
        dfc.append(read_file(file))
    df = pd.concat(dfc, axis=0, ignore_index=True)
    
    ipcLOG = np.log(df['Ipc'].astype('float64').replace(0.0,1.0))
    
    df['Ipc'] = ipcLOG
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
    # newkeys = list(set(newkeys))
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

clock = time()

randpic = sample([x for x in range(1, 202)], 3)
# randpic = [1, 2, 3, 4, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 82, 83, 84, 85, 123, 124, 126, 164, 165]

path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\pepprtmt\\desc\\'
files = [path + f'PFEC_PEPPRpep_{x}_desc.xlsx' for x in randpic]
df = df_compile(files)
#df.to_excel(path+'peppr_tmt_pep_MAX_desc.xlsx', index=False)
print(df.shape)

#df = pd.read_excel('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\desc\\'+
#                   'peppr_tmt_pep_55_desc.xlsx')

m = Chem.MolFromSmiles('C')
vals = Descriptors.CalcMolDescriptors(m)
keycards = list(vals.keys())
#keycards = ['PC2', 'PC3']

#newkeys= keycards

descrs = pd.DataFrame()
descrs['keys'] = keycards

#key_acc = [scale(df, keycard) for keycard in keycards]
#descrs['total acc'] = [key[0] for key in key_acc]
#descrs['target acc'] = [key[1] for key in key_acc]
#descrs['acc rank'] = descrs['target acc'] + descrs['total acc']
#descrs.to_csv('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\key_rank.csv',index=False)
#descrs = pd.read_csv('\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\key_rank.csv')

#newkeys = descrs3['keys'].values

record = 0
acc_tot,acc_targ=0,0
lists = []
# scores = [164.2, 155.1, 153.0, 160.4, 153.8, 153.8, 160.4, 160.0, 162.1, 160.3, 162.3, 164.2, 156.7, 160.6, 164.9, 162.1, 163.3, 164.9, 166.0, 163.5, 163.1, 164.4, 161.6, 165.2]
scores = [164.2, 160.4, 160.4, 160.0, 162.1, 160.3, 162.3, 164.2, 160.6, 164.9, 162.1, 163.3, 164.9, 166.0, 163.5, 163.1, 164.4, 161.6, 165.2]



# newkeys = ['FpDensityMorgan3', 'Chi2n', 'Chi4v', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'MaxAbsEStateIndex']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(4)
# lists.append(newkeys)

# # newkeys = ['MinEStateIndex', 'EState_VSA4', 'MaxEStateIndex']
# # # acc_tot, acc_targ = SVM_test(df, newkeys)
# # record = max(record, acc_tot+acc_targ)
# # print(6)
# # lists.append(newkeys)

# # newkeys = ['fr_diazo', 'EState_VSA4', 'MaxEStateIndex']
# # # acc_tot, acc_targ = SVM_test(df, newkeys)
# # record = max(record, acc_tot+acc_targ)
# # print(7)
# # lists.append(newkeys)

# newkeys = ['EState_VSA4', 'MaxAbsEStateIndex', 'VSA_EState9', 'FpDensityMorgan1', 'FpDensityMorgan3', 'BCUT2D_MRHI', 'AvgIpc', 'BCUT2D_MWLOW']#, 'MinEStateIndex', 'fr_alkyl_halide', 'fr_Ndealkylation2', 'NumSaturatedRings']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(8)
# lists.append(newkeys)

# # newkeys = ['FpDensityMorgan1', 'BCUT2D_MRHI', 'EState_VSA4', 'AvgIpc']
# # # acc_tot, acc_targ = SVM_test(df, newkeys)
# # record = max(record, acc_tot+acc_targ)
# # print(9)
# # lists.append(newkeys)


# # newkeys = ['FpDensityMorgan1', 'BCUT2D_MRHI', 'EState_VSA4', 'AvgIpc']
# # # acc_tot, acc_targ = SVM_test(df, newkeys)
# # record = max(record, acc_tot+acc_targ)
# # print(11)
# # lists.append(newkeys)


# newkeys = ['EState_VSA4', 'MaxAbsEStateIndex', 'VSA_EState9', 'FpDensityMorgan1', 'FpDensityMorgan3', 'BCUT2D_MRHI', 'AvgIpc', 'BCUT2D_MWLOW']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(13)
# lists.append(newkeys)
# #
# newkeys = ['EState_VSA4', 'MaxAbsEStateIndex', 'VSA_EState9', 'FpDensityMorgan1', 'FpDensityMorgan3', 'BCUT2D_MRHI', 'AvgIpc', 'MinEStateIndex']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(14)
# lists.append(newkeys)

# newkeys = ['fr_amidine', 'MaxEStateIndex', 'BalabanJ', 'EState_VSA4', 'fr_unbrch_alkane', 'FpDensityMorgan3']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(15)
# lists.append(newkeys)

# newkeys = ['BCUT2D_MRHI', 'VSA_EState9', 'EState_VSA4', 'FpDensityMorgan1', 'FpDensityMorgan3', 'MaxAbsEStateIndex']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(16)
# lists.append(newkeys)

# newkeys = ['fr_para_hydroxylation', 'NumAliphaticRings', 'fr_amide', 'BalabanJ', 'FpDensityMorgan1', 'FpDensityMorgan3']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(18)
# lists.append(newkeys)


# newkeys = ['FpDensityMorgan3', 'Chi2n', 'Chi4v', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'MaxAbsEStateIndex']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(20)
# lists.append(newkeys)

# # newkeys = ['NumSaturatedRings', 'FpDensityMorgan3', 'AvgIpc', 'fr_alkyl_halide', 'fr_Ndealkylation2']
# # # acc_tot, acc_targ = SVM_test(df, newkeys)
# # record = max(record, acc_tot+acc_targ)
# # print(21)
# # lists.append(newkeys)

# newkeys = ['SMR_VSA3', 'FpDensityMorgan1', 'MaxAbsEStateIndex', 'fr_lactone', 'FpDensityMorgan2']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(22)
# lists.append(newkeys)

# newkeys = ['MaxAbsEStateIndex', 'PEOE_VSA10', 'fr_Imine']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(23)
# lists.append(newkeys)

# newkeys = ['PEOE_VSA10', 'AvgIpc', 'fr_priamide', 'fr_thiocyan']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(24)
# lists.append(newkeys)

# newkeys = ['EState_VSA11', 'PEOE_VSA10']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(25)
# lists.append(newkeys)

# newkeys = ['MaxAbsEStateIndex', 'PEOE_VSA10']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(26)
# lists.append(newkeys)

# newkeys = ['PEOE_VSA10', 'PEOE_VSA13', 'fr_para_hydroxylation', 'MaxAbsEStateIndex']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(27)
# lists.append(newkeys)

# newkeys = ['FpDensityMorgan1', 'MaxEStateIndex', 'fr_unbrch_alkane', 'PEOE_VSA10', 'fr_C_S']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(28)
# lists.append(newkeys)


# newkeys = ['FpDensityMorgan3', 'PEOE_VSA10', 'fr_piperdine', 'fr_unbrch_alkane', 'fr_guanido', 'BCUT2D_MWLOW']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(29)
# lists.append(newkeys)

# newkeys = ['FpDensityMorgan3', 'Chi2v', 'MaxAbsEStateIndex', 'Chi4v', 'fr_term_acetylene']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(30)
# lists.append(newkeys)

# newkeys = ['MaxEStateIndex', 'fr_para_hydroxylation', 'AvgIpc', 'Chi3v', 'fr_ester']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(31)
# lists.append(newkeys)

# newkeys = ['Chi2v', 'fr_oxazole', 'FpDensityMorgan2', 'fr_unbrch_alkane']
# # acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(32)
# lists.append(newkeys)


# newkeys = ['BCUT2D_MWLOW', 'BCUT2D_MRHI', 'FpDensityMorgan1', 'BCUT2D_LOGPHI', 'fr_Nhpyrrole', 'BCUT2D_CHGLO']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(3)


# newkeys = ['MinEStateIndex', 'BCUT2D_MRHI', 'BCUT2D_MWLOW', 'FpDensityMorgan2', 'Chi0n']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(5)


# newkeys = ['BCUT2D_MRHI', 'MinEStateIndex', 'BCUT2D_MWLOW', 'fr_amide', 'FpDensityMorgan2', 'BCUT2D_MRLOW', 'fr_C_O_noCOO', 'FractionCSP3']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(6)


# newkeys = ['MinEStateIndex', 'BCUT2D_CHGLO', 'Chi4v', 'fr_guanido', 'AvgIpc', 'fr_ether', 'fr_unbrch_alkane']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(8)


# newkeys = ['MaxEStateIndex', 'AvgIpc', 'FpDensityMorgan3', 'MinEStateIndex', 'PEOE_VSA4', 'fr_C_O_noCOO', 'fr_azo']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(11)


newkeys = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']
acc_tot, acc_targ = SVM_test(df, newkeys)
record = max(record, acc_tot+acc_targ)
print(3)

# newkeys = ['fr_para_hydroxylation', 'qed', 'fr_C_O_noCOO', 'fr_unbrch_alkane', 'fr_priamide', 'FpDensityMorgan3']
# acc_tot, acc_targ = SVM_test(df, newkeys)
# record = max(record, acc_tot+acc_targ)
# print(3)

# for _ in range(20):
#     keys = sample(newkeys,sample([17,16,15,14,13, 12, 11, 10, 9, 8, 7, 6, 5], 1)[0])
#     acc_tot, acc_targ = SVM_test(df, keys)
#     record = max(record, acc_tot+acc_targ)
#     print('-')

#newkeys = scale_maxer(df, descrs, [descrs.sort_values('acc rank')['keys'].values[-1]], 8, 3)
#newkeys = scale_maxer(df, descrs, newkeys, 8, 3)



print()
print(record)
print(f'Total runtime: {round(time()-clock, 1)} sec')
