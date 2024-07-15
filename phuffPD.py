# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:43:47 2022

@author: HUFFMP
"""
print('Loading Python Modules: ', end='')
import pandas as pd
print('', end='-')
import math
print('', end='-')
import numpy as np
print('', end='-')
import matplotlib.pyplot as plt
print('', end='-')
import os
print('', end='-')
from sklearn.decomposition import PCA
print('', end='-')
import matplotlib.patches as mpatches
print('', end='-')
import sys
print('', end='-')
import pathlib
print('', end='-')
from time import sleep
print('', end='-')
pd.options.mode.chained_assignment = None
print('-')

def crash():
    print('\nExiting phuffPD. \n\n')
    sleep(1.2)
    sys.exit()


def tmtPD(df, sav):
    # newpath = f'\\\\amer.pfizer.com\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\Library\\external\\{sav}'
    # if not os.path.exists(newpath):
    #     os.makedirs(newpath)
    newpath = sav

    ab_target = [x for x in df.columns if str(x).count('Abundance:')>0]

    tmtVals = df[ab_target].dropna()
    tmtVals['TEMP'] = tmtVals.mean(axis=1)

    q15 = tmtVals['TEMP'].quantile(0.15)
    q85 = tmtVals['TEMP'].quantile(0.85)

    tmtVals = tmtVals[tmtVals['TEMP']>q15]
    tmtVals = tmtVals[tmtVals['TEMP']<q85]
    tmtVals = tmtVals.reset_index(drop=True)
    del tmtVals['TEMP']


    plt.boxplot(tmtVals, showfliers=True)
    plt.ylabel('Abundance')
    plt.xlabel('TMT Lane')
    plt.title('TMT Lanes: Boxplot with Fliers')
    plt.xticks([i-1 for i in range(1, 19)], [str(i) for i in range(1, 19)])
    plt.savefig(newpath+'/box.png', dpi=1200)
    print(f'boxplot saved to folder')
    plt.clf()

    plt.boxplot(tmtVals, showfliers=False)
    plt.ylabel('Abundance')
    plt.xlabel('TMT Lane')
    plt.title('TMT Lanes: Boxplot')
    plt.xticks([i-1 for i in range(1, 19)], [str(i) for i in range(1, 19)])
    plt.savefig(newpath+'/box_simple.png', dpi=1200)
    print(f'simple boxplot saved to folder')
    plt.clf()

    plt.bar(tmtVals.columns, tmtVals.mean(axis=0).values)
    plt.ylabel('Mean Abundance')
    plt.xlabel('TMT Lane')
    plt.xticks([i-1 for i in range(1, 19)], [str(i) for i in range(1, 19)])
    plt.title('TMT Lanes: Mean Abundances')
    plt.savefig(newpath+'/bar.png', dpi=1200)
    print(f'barplot saved to folder')
    plt.clf()

    tmtMXR = pd.DataFrame(index=range(1, 19))
    tmtMXR['means'] = tmtVals.mean().values

    means = tmtVals.mean()
    meanmin = max(means.values)
    for val in means.values:
        if (val < meanmin) and (val*25 > max(means.values)):
            meanmin = val

    tmtMXR['rel. ratio'] = tmtMXR['means'] / meanmin
    tmtMXR['add. percent'] = 1 / tmtMXR['rel. ratio']

    tmtMXR['add. percent'][tmtMXR['add. percent']>1] = 0
    tmtMXR = tmtMXR.round({'add. percent':3, 'rel. ratio':4, 'means':4})

    tmtMXR.columns = ['Abundance Means', 'Relative Abundance Ratio', 'New Mix Ratio']

    plt.table(cellText=tmtMXR.values, colLabels=tmtMXR.columns, loc='center', cellLoc='center')
    plt.axis('off')
    plt.title('TMT mix ratio calculation', y=1.08)
    plt.tight_layout()
    plt.xticks([])
    plt.yticks([])
    plt.savefig(newpath+'/TMTsummary.png', dpi=1200)
    tmtMXR.to_excel(newpath+'/TMTsummary.xlsx')
    print(f'TMT mix ratio table saved to folder')
    plt.clf()
    return ab_target


def ptna_name_rep(name):
    if name.count('_') > 0:
        if name[-1:] != '_':
            return ptna_name_rep(name[:-1])
    else:
        name += '_'
    return name[:-1]


def ptna_ellipsewidth(pointrem, pts):
    m = (pts[1][1]-pts[0][1]) / (pts[1][0]-pts[0][0])
    dx = abs((pointrem[1]) - (pts[0][1]+((pointrem[0]-pts[0][0]) * m)))
    tridist = dx * math.sin(90-math.atan(m))
    return tridist


def pcaPD(df, sav, numreps):
    # newpath = f'\\\\amer.pfizer.com\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\Library\\external\\{sav}'
    # if not os.path.exists(newpath):
    #     os.makedirs(newpath)
    newpath = sav

    abundances = [x for x in df.columns if str(x).count('Abundance:')>0]

    if numreps == 3:
        names = ['1', '1', '1', '2', '2', '2',
                 '3', '3', '3', '4', '4', '4',
                 '5', '5', '5', '6', '6', '6']
    elif numreps == 2:
        names = ['1', '1', '2', '2', '3', '3',
                 '4', '4', '5', '5', '6', '6',
                 '7', '7', '8', '8', '9', '9']
    else:
        names = ['1', '2', '3', '4', '5', '6',
                 '7', '8', '9', '10', '11', '12',
                 '13', '14', '15', '16', '17', '18']

#    print(numreps)


    xf = df[df['# PSMs']>4][abundances].dropna().replace(0.0, 1.0)

    xf['ab sum'] = df[abundances].dropna().sum(axis=1)
    xf = xf[xf['ab sum']>xf['ab sum'].quantile(0.2)]
    del xf['ab sum']

    df_main_transformed = np.log2(xf)
    
    df_col_median = df_main_transformed.median()
    df = df_main_transformed-df_col_median

    ncolors = []
    colorchoices = ['blue', 'orange', 'green', 'red', 'purple',
                    'brown', 'pink', 'gray', 'olive', 'cyan',
                    'navy', 'lightsalmon', 'lime', 'mediumturquoise', 'maroon',
                    'skyblue', 'fuchsia', 'darkslateblue']

    color_count = 0
    color_dict = {}

    for name in names:
        if ptna_name_rep(name) in color_dict.keys():
            ncolors.append(color_dict[ptna_name_rep(name)])
        else:
            color_dict[ptna_name_rep(name)] = colorchoices[color_count]
            ncolors.append(colorchoices[color_count])
            color_count += 1


    x = df.values
#        df.to_csv('pca_input.csv', index=False)
    pca = PCA(n_components=2)

    pCones = pca.fit_transform(x.T)
    principalDf = pd.DataFrame(data=pCones, columns = ['PC1', 'PC2'])

    xs = principalDf['PC1']
    ys = principalDf['PC2']
    # print(xs)

    ncolors = pd.Series(ncolors)

    plt.clf()
    plt.close()

    plt.scatter(xs, ys, c=np.array(ncolors))

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')

    plt.legend(handles=[mpatches.Patch(color=col, label=text) for text, col in color_dict.items()],
                        prop={'size': 6})
    plt.title('PCA of Lane Similarity')
    plt.savefig(newpath+f'/TMT_pca{numreps}.png', dpi=1200)
    plt.clf()
    plt.close()


switch=True
while switch:
    print('Welcome to phuffPD. Proteomics Summaries on the fly.')
    filename = input('Filename? : ')
    fileSuffix = pathlib.Path(filename).suffix
    foldersav = pathlib.Path(filename).stem
    print('Searching for and reading file.')
    if fileSuffix == '.xlsx':
        df = pd.read_excel(filename)
        switch=False
    elif fileSuffix == '.csv':
        df = pd.read_csv(filename)
        switch=False
        print('File read successfully.')
    else:
        print('File type not recognized.')
        switch = True
    print('File read successfully.')


sav = filename[:-(len(fileSuffix)+len(foldersav))]
newpath = f'{sav}\\PHpca_{foldersav}'
if not os.path.exists(newpath):
    os.makedirs(newpath)

df = df.dropna(thresh=2)

abund = tmtPD(df, newpath)


pcaPD(df, newpath, 3)
pcaPD(df, newpath, 2)
pcaPD(df, newpath, 1)



crash()