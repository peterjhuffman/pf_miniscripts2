# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:20:28 2023

@author: HUFFMP
"""

# Imports
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn import svm
from sklearn.metrics import confusion_matrix
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.decomposition import PCA
from time import time
from random import sample
import pickle
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier

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
    nclock = timer_quant(time()-timesplit)
    timesplit = time()
    print(message.ljust(55)+f"{nclock[0]} hr {nclock[1]} min {nclock[2]} sec")
    return timesplit


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
    return df


def df_compile(list_files):
    cols = ['Modifications','Positions in Master Proteins','Sequence','Accession','Acquired','SMILES string',
            'fr_Nhpyrrole','fr_guanido','BCUT2D_CHGLO','MinEStateIndex','Chi4v','fr_ether','FpDensityMorgan1',
            'AvgIpc','fr_unbrch_alkane','BCUT2D_MRHI','BCUT2D_MWLOW','BCUT2D_LOGPHI']
    dfc = []
    for file in list_files:
        print(file.split('\\')[-1].split('.')[0][:-5], end='}  ')
        dfc.append(read_file(file)[cols])
    df = pd.concat(dfc, axis=0, ignore_index=True, sort=True)
    
    if 'Ipc' in df.columns:
        ipcLOG = np.log(df['Ipc'].astype('float64').replace(0.0,1.0))
        df['Ipc'] = ipcLOG


    return df.dropna(subset=cols[6:])

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


def SVM_multifile(df1, df2, newkeys):
    print(str(len(newkeys))+' descriptors analyzed')
    
    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))
    
    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto', probability=True)
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

    return clf


def SVM_monofile(df1, newkeys):
    print(str(len(newkeys))+' descriptors analyzed')
    
#    # Generate data
#    X_train, X_test, y_train, y_test = train_test_split(df1[newkeys].values, df1['Acquired'].values, test_size=0.25)
#    bal = 1/(y_train.sum()/len(y_train))
#    
#    # Initialize SVM classifier
#    wetweight=1.0
#    weights = {0:1.0, 1:bal*wetweight}
#    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto', probability=True)
#    clf = clf.fit(X_train, y_train)
#    
#    
#    # Predict the test set
#    predictions = clf.predict(X_test)
#    
#    
#    # Generate confusion matrix
#    matrix = confusion_matrix(y_test, predictions, labels=clf.classes_)
#    
#    print(matrix)
#    acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
#    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)
#    
#    print(f'Total Accuracy: {acc_tot}%')
#    print(f'Target Accuracy: {acc_targ}%')
#    print(f'Accuracy Rank: {acc_targ+acc_tot}')


    bal = 1/(df1['Acquired'].values.sum()/len(df1['Acquired'].values))
    
    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel='rbf', class_weight=weights, gamma='auto', probability=True)
    clf2 = clf.fit(df1[newkeys].values, df1['Acquired'].values)

    return clf2


def lin_SVM_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = svm.SVC(kernel="linear", C=0.025, random_state=42, class_weight=weights)
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

    return clf


def gaus_proc_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = GaussianProcessClassifier(1.0 * RBF(1.0), random_state=42)
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

    return clf


def dec_tree_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = DecisionTreeClassifier(max_depth=5, random_state=42, class_weight=weights)
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

    return clf


def rand_forest_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf =  RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1, random_state=42, class_weight=weights)
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

    return clf


def mlp_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf =  MLPClassifier(alpha=1, max_iter=1000, random_state=42)
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

    return clf


def ada_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = AdaBoostClassifier(random_state=42)
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

    return clf


def nativ_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
    weights = {0:1.0, 1:bal*wetweight}
    clf = GaussianNB()
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

    return clf


def knn_multifile(df1, df2, newkeys):
    df1.to_csv('df1.csv')
    print(newkeys)

    print(str(len(newkeys))+' descriptors analyzed')

    # Generate data
    X_train, X_test, y_train, y_test = df1[newkeys].values, df2[newkeys].values, df1['Acquired'].values, df2['Acquired'].values
    bal = 1/(y_train.sum()/len(y_train))

    # Initialize SVM classifier
    wetweight=1.0
#    weights = {0:1.0, 1:bal*wetweight}
    clf = KNeighborsClassifier(3)
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

    return clf


def random_key(keys, sample_v):
    return sample([x for x in keys], sample_v)


def clf_test(df, clf, newkeys):
    x, y = df[newkeys].values, df['Acquired'].values
    predictions = clf.predict(x)
    df['predictions'] = predictions
#    df['predscore'] = pd.DataFrame(clf.predict_proba(x))[1][:]
    df.to_csv('predictions on HELA.csv', index=False)
    matrix = confusion_matrix(y, predictions, labels=clf.classes_)
    
    print(matrix)
    acc_tot = round(((matrix[0][0]/matrix[0].sum())+(matrix[1][1]/matrix[1].sum()))*50, 1)
    acc_targ = round(matrix[1][1]/max(1,matrix[1].sum())*100, 1)
    
    print(f'Total Accuracy: {acc_tot}%')
    print(f'Target Accuracy: {acc_targ}%')
    print(f'Accuracy Rank: {acc_targ+acc_tot}')

    return acc_tot, acc_targ

clock = time()
globkeys = ['fr_Nhpyrrole', 'fr_guanido', 'BCUT2D_CHGLO', 'MinEStateIndex', 'Chi4v', 'fr_ether', 'AvgIpc', 'fr_unbrch_alkane', 'BCUT2D_MWLOW', 'BCUT2D_MRHI', 'BCUT2D_LOGPHI', 'FpDensityMorgan1']

m = Chem.MolFromSmiles('C')
vals = Descriptors.CalcMolDescriptors(m)
keycards = list(vals.keys())
#keycards = ['PC2', 'PC3']

#newkeys= keycards
timesplit = time()
descrs = pd.DataFrame()
descrs['keys'] = keycards



#r1 = sample([x for x in range(0, 201)], 10)#200
r1 = sample([x for x in range(0,78)], 10)
path1 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\chymo TRPM8\\bESI training\\desc\\'
files1 = [path1 + f'PFAS_CHYpep_{x}_desc.xlsx' for x in r1]

#r2 = sample([x for x in range(0, 201)], 10)
#path2 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\pepprtmt\\desc\\'
#files2 = [path2 + f'PFEC_PEPPRpep_{x}_desc.xlsx' for x in r2]
#r2 = sample([x for x in range(0, 33)], 10)#31
##r2 = [5]
#path2 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI TMTyl-test\\desc\\'
#files2 = [path2 + f'PFAS_MASTpep_{x}_desc.xlsx' for x in r2]
##files2 = ['\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\examplepeps\\desc\\ar5deg_sample_desc.xlsx']
#
#
#r3 = sample([x for x in range(0, 94)], 9)#92
#path3 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI ascendhela-test\\desc\\'
#files3 = [path3 + f'PFEC_ascHELApep_{x}_desc.xlsx' for x in r3]
##files1 += files3
#
#r4 = sample([x for x in range(0, 89)], 9)#82
#if 9 in r4:
#    r4.remove(9)
#path4 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI hela-test\\desc\\'
#files4 = [path4 + f'PFEC_HELApep_{x}_desc.xlsx' for x in r4]
#files1 += files4
#
#r5 = sample([x for x in range(0, 14)]+[x for x in range(100, 102)]+[x for x in range(200, 217)]+[x for x in range(300, 305)]+
#             [x for x in range(400, 407)], 14)#30
#path5 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI tmtrtsRG-test\\desc\\'
#files5 = [path5 + f'PFAS_RTSpep_{x}_desc.xlsx' for x in r5]
#files1 += files5
#
#r6 = sample([x for x in range(210, 243)], 9)#31
#path6 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI TMTyl-test\\desc\\'
#files6 = [path6 + f'PFAS_MASTpep_{x}_desc.xlsx' for x in r6]
#files1 += files6
#files1 += files2
#
#r7 = sample([x for x in range(0, 15)], 9)#14
#path7 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI TMT2-test\\desc\\'
#files7 = [path7 + f'PFAS_E3Opep_{x}_desc.xlsx' for x in r7]
#files1 += files7
#
#r8 = sample([x for x in range(0, 34)], 9)#33
#files8 = [path7 + f'PFAS_YLpep_{x}_desc.xlsx' for x in r8]
#files1 += files8
#
#r9 = sample([x for x in range(0, 63)], 9)#75
#path9 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\bESI TMT-test\\desc\\'
#files9 = [path9 + f'PFAS_LFQpep_{x}_desc.xlsx' for x in r9]
#files1 += files9
#
#r0 = sample([x for x in range(30, 76)], 9)#34
#files0 = [path9 + f'PFAS_TMTpep_{x}_desc.xlsx' for x in r0]
#files1 += files0
#
#r0 = sample([x for x in range(80, 183)], 9)#101
#files0 = [path9 + f'PFAS_XMPpep_{x}_desc.xlsx' for x in r0]
#files1 += files0


#files1.reverse()
#df2 = df_compile(files2)
#print(df2.shape)
print()
print()
#df1 = df_compile(files1)
df1 = df_compile(files1)
print(df1.shape)




#
#print('\n\nFiles Loaded. Commencing Training.\n')
timesplit = timed_print('Input loaded.', timesplit)
#
#clf = SVM_monofile(df1, globkeys)
#with open(path1+'bESI_v3_2_chymomax.pkl', 'wb') as fid:
#    pickle.dump(clf, fid)
#timesplit = timed_print('SVM chymo', timesplit)
#print()
#
#clf = rand_forest_multifile(df1, df2, globkeys)
#timesplit = timed_print('rand forest', timesplit)
#print()
###
#clf = dec_tree_multifile(df1, df2, globkeys)
#timesplit = timed_print('dec tree', timesplit)
#print()
#
#clf = lin_SVM_multifile(df1, df2, globkeys)
#timesplit = timed_print('lin svm', timesplit)
#print()
#
#clf = SVM_multifile(df1, df2, globkeys)
#timesplit = timed_print('svm', timesplit)
#print()
#
#
#
#
#with open(path1[:-14]+'bESI_v2_6_0p.pkl', 'wb') as fid:
#    pickle.dump(clf, fid)

v1 = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\bESI\\'
#
#paths = ['bESI_v2_0.pkl', 'bESI_v2_2.pkl', 'bESI_v2_3_1.pkl', 'bESI_v2_4_0p.pkl', 'bESI_v2_4_1.pkl', 'bESI_v2_4_2.pkl',
#         'bESI_v2_4_3p.pkl', 'bESI_v2_4_decTree.pkl', 'bESI_v2_4_linSVM.pkl', 'bESI_v2_4_rforest.pkl', 'bESI_v2_4_svm.pkl',
#         'bESI_v2_5_2dt.pkl','bESI_v2_5_2rf.pkl','bESI_v2_5_rf.pkl']
#
paths = ['bESI_v3_2_chymomax.pkl']
#paths = ['bESI_v2_6_0.pkl']
for pathy in paths:
    timesplit = timed_print(pathy, timesplit)
    with open(v1+pathy, 'rb') as fid:
        clf_R = pickle.load(fid)
    clf_test(df1, clf_R, globkeys)
    print()




#print()
print(f'Total runtime: {round(time()-clock, 1)} sec')
