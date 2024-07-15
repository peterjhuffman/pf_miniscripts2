# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:30:39 2024

@author: HUFFMP
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


print('Searching for file...')
path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\python scripts\\sizec ptmx\\'
#filename = 'nSEC_atlas_v1_normalized.xlsx'
filename = 'nSEC_atlas_v1_cleaned.xlsx'

df = pd.read_excel(path+filename)
print('File read.')

# -----------------------------------------------------------------------------------------
METHOD = 'TARGET'
SAVE = False
# -----------------------------------------------------------------------------------------
#targets = ['CDK1', 'CDK11B', 'CDK12', 'CDK13', 'CDK16', 'CDK17','CDK18','CDK2','CDK4','CDK5','CDK6',
#           'CDK7','CDK9','CDKAL1','CDKN1B','CDKN2A','CCNA2','CCNB1','CCNB2','CCND1','CCND3','CCNH','CCNK','CCNL1','CCNT1','CCNY']
#targets = ['GALE', 'TXN','GPX4','RASA2','NOLC1','GMPPB','PLCB4','ADAP1',
#           'S100A14','BAG5','GMPPA','C1orf52','SNX15','GBA','C15orf40',
#           'PADI1','BBOX1','WIZ','SERPINB9','RPS23','SLC25A3','RPS4X',
#           'SLC25A5','RPS11','RPS8','SLC25A6','ATP5J2','XPO1','ATP5C1'
#           ,'RPS16','RPL6','RPS18','RPS13','TMX1']
targets= ['ZNF326']
phmap = 'viridis'


if METHOD == 'TARGET':
    for tg in targets:
        dfs = pd.DataFrame({'Model':df['Model'], 'Fraction':df['Fraction'], 'Target':df[tg]})
        dfs = dfs[dfs['Fraction']!=1].reset_index(drop=True)
        dfs = dfs[dfs['Fraction']!=24].reset_index(drop=True)

        dfscells = [x for x in dfs['Model'].unique()]
        dictresh = {'Fraction':dfs[dfs['Model']==dfscells[0]].reset_index(drop=True)['Fraction']}

        for cell in dfscells:
            dictresh[cell]=dfs[dfs['Model']==cell].reset_index(drop=True)['Target']

        dfr = pd.DataFrame(dictresh)

        print(dfr.shape)
        print(dfr.columns)

        plt.figure(figsize=(20,12))
        for cell in dfscells:
            plt.plot(dfr['Fraction'], dfr[cell], label=cell)
        plt.legend()
        plt.title(f'{tg} sizec')
        plt.xlabel('Fraction #')
        plt.ylabel('Abundance')
        if SAVE:
            plt.savefig(path+f'pulls\\{tg}_line.png',dpi=200)
#        plt.show()
        plt.clf()

        plt.figure(figsize=(15,9))
        plt.imshow(dfr[dfscells],cmap=phmap)
        cbar= plt.colorbar()
        cbar.set_label('Abundance')
        plt.title(f'{tg} sizec')
        plt.ylabel('Fraction #')
        plt.xlabel('Cell Line')
        plt.yticks(ticks=[-0.5, 0, 5, 10, 15, 21, 21.5], labels=['', 2, 7, 12, 17, 23, ''])
        plt.xticks([i-1 for i in range(1, len(dfscells)+1)], dfscells, rotation=90)
        if SAVE:
            plt.savefig(path+f'pulls\\{tg}_heat.png',dpi=200, bbox_inches='tight')
#        plt.show()
        plt.clf()
#
        newlabels = [f"{cell} cellnorm" for cell in dfscells]
        dfrn = dfr.copy()
        for cell in dfscells:
            cellstat = sum(dfrn[cell])
            dfrn[f"{cell} cellnorm"] = dfrn[cell]/max(cellstat, 0.1)
        
        plt.figure(figsize=(15,9))
        plt.imshow(dfrn[newlabels],cmap=phmap)
        cbar=plt.colorbar()
        cbar.set_label('Abundance (norm)')
        plt.title(f'{tg} sizec')
        plt.ylabel('Fraction #')
        plt.xlabel('Cell Line')
        plt.yticks(ticks=[-0.5, 0, 5, 10, 15, 21, 21.5], labels=['', 2, 7, 12, 17, 23, ''])
        plt.xticks([i-1 for i in range(1, len(dfscells)+1)], dfscells, rotation=90)
#        if SAVE:
#            plt.savefig(path+f'pulls\\{tg}_heatn.png',dpi=200, bbox_inches='tight')
#        plt.show()
        plt.clf()

        x = dfrn[newlabels].values
        pca = PCA(n_components=2)
        pCones = pca.fit_transform(x.T)
        subDf = pd.DataFrame(data=x, columns=newlabels)
        principalDf = pd.DataFrame(data=pCones, columns = ['PC1', 'PC2'])
        principalDf['Cell']=dfscells
        xs = principalDf['PC1']
        ys = principalDf['PC2']
        plt.scatter(xs, ys)
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('PCA of Lane Similarity')
        if SAVE:
            plt.savefig(path+f'pulls\\{tg}_pca.png',dpi=200)
#        plt.show()
        plt.clf()


        x = dfrn[newlabels].values
        pca = PCA(n_components=1)
        pCones = pca.fit_transform(x.T)

        dfpc = dfrn[newlabels]
        dfpc.loc[22, :] = [x[0] for x in pCones]
        dfpca = dfpc.sort_values(by=22, axis=1)
        dfpca=dfpca.drop(22)
        
        plt.figure(figsize=(15,9))
        plt.imshow(dfpca,cmap=phmap)
        cbar=plt.colorbar()
        cbar.set_label('Abundance (norm)')
        plt.title(f'{tg} sizec')
        plt.ylabel('Fraction #')
        plt.xlabel('Cell Line')
        plt.yticks(ticks=[-0.5, 0, 5, 10, 15, 21, 21.5], labels=['', 2, 7, 12, 17, 23, ''])
        plt.xticks([i-1 for i in range(1, len(dfscells)+1)], [x[:-9] for x in dfpca.columns], rotation=90)
        if SAVE:
            plt.savefig(path+f'pulls\\{tg}_heatar.png',dpi=200, bbox_inches='tight')
        plt.show()
        plt.clf()
























































