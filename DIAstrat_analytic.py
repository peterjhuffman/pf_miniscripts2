# -*- coding: utf-8 -*-
"""
Created on Wed May  8 09:01:50 2024

@author: HUFFMP
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from sklearn import datasets, linear_model
from sklearn.metrics import r2_score
from sklearn.decomposition import PCA
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from sklearn.impute import KNNImputer

def knn_tn_impute(X, k=5):
    """
    Impute missing values using KNN-TN.
    
    Args:
        X (numpy.ndarray): Input data matrix with missing values.
        k (int): Number of neighbors to consider (default: 5).
    
    Returns:
        numpy.ndarray: Imputed data matrix.
    """
    imputer = KNNImputer(n_neighbors=k)
    imputed_X = imputer.fit_transform(X)
    return imputed_X



def colorsplit(pct):
    rgb_color = (pct, 0, (1-pct))
    hex_color = mcolors.to_hex(rgb_color)
    return str(hex_color)


def ptna_name_rep(name):
    if name.count('_') > 0:
        if name[-1:] != '_':
            return ptna_name_rep(name[:-1])
    else:
        name += '_'
    return name[:-1]


def phuff_impute(rs, ls, pcam, exvar):

    nanloc = rs[rs.isna()].index[0]
    truloc = [x for x in rs[~rs.isna()].index]

    val1 = rs[truloc[0]]
    val2 = rs[truloc[1]]

    dist1 = [(x-y)**2 for x,y in zip(pcam.loc[nanloc], pcam.loc[truloc[0]])]
    dist2 = [(x-y)**2 for x,y in zip(pcam.loc[nanloc], pcam.loc[truloc[1]])]

    var1 = sum([x*y for x,y in zip(dist1, [i for i in exvar])])
    var2 = sum([x*y for x,y in zip(dist2, [i for i in exvar])])

    rat1 = var1 / (var1+var2)
    rat2 = var2 / (var1+var2)

    rs.loc[nanloc] = (rat1*val1)+(rat2*val2)

    return rs


path = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\34. DIA practice\\exports\\'

print('Reading Files.')

# prot = pd.read_excel(path+'report.pg_matrix_1miscut.xlsx')
# pep = pd.read_excel(path+'report.pr_matrix_1mixcut.xlsx')

# prot.columns = [x for x in prot.columns[:5]] + [f"LFQ {x+1}" for x in range(20)]
# pep.columns = [x for x in pep.columns[:10]] + [f"LFQ {x+1}" for x in range(20)]

# protval = prot[prot.columns[5:]]

print('Files Read.')


# MONTE CARLO ADDITIVE INTERSECTION =======================================================================================

# sd = {}
# for x in range(len(protval.columns)):
#     s = x+1
#     l = [] 
#     for rep in range(40):
#         rs = random.sample([x for x in protval.columns], s)
#         l.append(protval.dropna(subset=rs).shape[0])
#     sd[s] = np.median(l)
#     print('',end='-')

# np.mean([len(protval[x].dropna()) for x in protval.columns])
# plt.plot(sd.keys(), sd.values())
# # plt.ylim(0,5700)
# plt.xlabel('replicates')
# plt.ylabel('shared protein IDs')
# plt.title('shared prot IDs by #rep')
# plt.savefig(path+'depict\\protINTERSECT_01.png',dpi=300)
# plt.show()
# plt.clf()

# plt.plot(sd.keys(), sd.values())
# plt.ylim(0,5700)
# plt.xlabel('replicates')
# plt.ylabel('shared protein IDs')
# plt.title('shared prot IDs by #rep')
# plt.savefig(path+'depict\\protINTERSECT_02.png',dpi=300)
# plt.show()
# plt.clf()


# r2 of proteome mix =============================================================================================

# a = []
# b = []

# for i in protval.index:
#     protslice = protval.iloc[i].dropna()
#     if len(protslice)>1:
#         for _ in range(len(protslice)*2):
#             rs = random.sample([x for x in protslice.values], 2)
#             a.append(rs[0])
#             b.append(rs[1])

# plt.scatter(a, b, s=0.08)


# x = np.asarray(a).reshape(-1,1)
# y = np.asarray(b).reshape(-1,1)

# regr = linear_model.LinearRegression()
# regr.fit(x, y)
# score = r2_score(y, regr.predict(x))
# print

# plt.ylabel('LFQ')
# plt.xlabel('LFQ')
# plt.title(f'LFQ replicability   -  r2: {round(score, 3)}')
# plt.savefig(path+'depict\\lfq_REPLICATE_01.png',dpi=300)
# plt.show()
# plt.clf()



# a = []
# b = []

# for i in protval.index:
#     protslice = protval.iloc[i].dropna()
#     if len(protslice)>1:
#         for _ in range(len(protslice)*2):
#             rs = random.sample([x for x in protslice.values], 2)
#             a.append(np.log2(rs[0]))
#             b.append(np.log2(rs[1]))

# plt.scatter(a, b, s=0.008)


# x = np.asarray(a).reshape(-1,1)
# y = np.asarray(b).reshape(-1,1)

# regr = linear_model.LinearRegression()
# regr.fit(x, y)
# score = r2_score(y, regr.predict(x))
# print

# plt.ylabel('LFQ log2')
# plt.xlabel('LFQ log2')
# plt.title(f'LFQ log2 replicability   -  r2: {round(score, 3)}')
# plt.savefig(path+'depict\\lfq_REPLICATE_02.png',dpi=500)
# plt.show()



# CV by LFQ ====================================================================================================


# a = []
# b = []

# cv = lambda x: np.std(x) / np.mean(x) * 100

# for i in protval.index:
#     protslice = protval.iloc[i].dropna()
#     cva = cv(protslice.values)
#     for v in protslice.values:
#         a.append(np.log2(v))
#         b.append(cva)
# # 
# plt.scatter(a, b, s=0.001)
# plt.ylabel('CV')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 replicability')
# plt.savefig(path+'depict\\lfq_CV_01.png',dpi=300)
# plt.show()

# c = []
# d = []
# d1 = []

# abdf = pd.DataFrame({'val':a, 'err':b})

# valmin = min(a)
# valmax = max(a)

# step = (valmax-valmin)/20


# while valmin<valmax:
#     ab1mask = abdf['val']>valmin
#     ab2mask = abdf['val']<valmin+step

#     c.append(np.mean([valmin, valmin+step]))
#     abdf_k = abdf[ab1mask&ab2mask]

#     d.append(abdf_k['err'].quantile(0.95))
#     d1.append(abdf_k['err'].median())
#     valmin+=step


# plt.plot(c, d1)
# plt.ylabel('CoV')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 median%CV')
# plt.savefig(path+'depict\\lfq_CVx_01.png',dpi=300)
# plt.show()

# plt.plot(c, d)
# plt.ylabel('CoV')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 95%CV')
# plt.savefig(path+'depict\\lfq_CVx_02.png',dpi=300)
# plt.show()




# %error by LFQ ====================================================================================================

# a = []
# b = []

# # protval['mean'] = protval.mean(axis=1)


# for i in protval.index:
#     protslice = protval.iloc[i].dropna()
#     mstat = np.mean(protslice.values)
#     for v in protslice.values:
#         a.append(np.log2(v))
#         b.append(100*abs(v-mstat)/v)

# plt.scatter(a, b, s=0.001)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 %error')
# plt.savefig(path+'depict\\lfq_mstat_01.png',dpi=300)
# plt.show()


# c = []
# d = []
# d1 = []


# abdf = pd.DataFrame({'val':a, 'err':b})

# valmin = min(a)
# valmax = max(a)

# step = (valmax-valmin)/100


# while valmin<valmax:
#     ab1mask = abdf['val']>valmin
#     ab2mask = abdf['val']<valmin+step

#     c.append(np.mean([valmin, valmin+step]))
#     abdf_k = abdf[ab1mask&ab2mask]

#     d.append(abdf_k['err'].quantile(0.95))
#     d1.append(abdf_k['err'].median())
#     valmin+=step

# plt.scatter(c, d1, s=1)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 median%error')
# plt.savefig(path+'depict\\lfq_mstat_02.png',dpi=300)
# plt.show()

# plt.scatter(c, d, s=1)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 95%error')
# plt.savefig(path+'depict\\lfq_mstat_03.png',dpi=300)
# plt.show()

# plt.plot(c, d1)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 median%error')
# plt.savefig(path+'depict\\lfq_mstat_04.png',dpi=300)
# plt.show()

# plt.plot(c, d)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 95%error')
# plt.savefig(path+'depict\\lfq_mstat_05.png',dpi=300)
# plt.show()

# plt.plot(c, d1)
# plt.ylim(0,100)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 median%error')
# plt.savefig(path+'depict\\lfq_mstat_06.png',dpi=300)
# plt.show()

# plt.plot(c, d)
# plt.ylim(0,100)
# plt.ylabel('%error')
# plt.xlabel('LFQ log2')
# plt.title('LFQ log2 95%error')
# plt.savefig(path+'depict\\lfq_mstat_07.png',dpi=300)
# plt.show()



# replicate ID boost ====================================================================================================

# sd = {}
# for x in range(len(protval.columns)):
#     s = x+1
#     l = []
#     for rep in range(40):
#         rs = random.sample([x for x in protval.columns], s)
#         l.append(protval.dropna(subset=rs, thresh=1).shape[0])
#     sd[s] = np.median(l)
#     print('',end='-')

# # np.mean([len(protval[x].dropna()) for x in protval.columns])
# plt.plot(sd.keys(), sd.values())
# # plt.ylim(0,5700)
# plt.xlabel('replicates')
# plt.ylabel('shared protein IDs')
# plt.title('shared prot IDs by #rep')
# plt.savefig(path+'depict\\protUNION_01.png',dpi=300)
# plt.show()
# plt.clf()

# plt.plot(sd.keys(), sd.values())
# plt.ylim(0,6000)
# plt.xlabel('replicates')
# plt.ylabel('shared protein IDs')
# plt.title('shared prot IDs by #rep')
# plt.savefig(path+'depict\\protUNION_02.png',dpi=300)
# plt.show()
# plt.clf()

# plt.plot(sd.keys(), sd.values())
# # plt.ylim(0,6000)
# plt.xlim(0.5,6.5)
# plt.xlabel('replicates')
# plt.ylabel('shared protein IDs')
# plt.title('shared prot IDs by #rep')
# plt.savefig(path+'depict\\protUNION_03.png',dpi=300)
# plt.show()
# plt.clf()



# avg missingval by lfq ====================================================================================================

# a = []
# b = []

# for i in protval.index:
#     protslice = protval.iloc[i].dropna()
#     mstat = np.mean(protslice.values)
#     a.append(np.log2(mstat))
#     b.append(len(protslice)*5)

# plt.scatter(a,b,s=0.1)
# plt.xlabel('LFQ log2')
# plt.ylabel('#missing vals')
# plt.title('#missing vals by lfq')
# plt.savefig(path+'depict\\protMISS_01.png',dpi=300)
# plt.show()
# plt.clf()

# c = []
# d = []
# d1 = []


# abdf = pd.DataFrame({'val':a, 'err':b})

# valmin = min(a)
# valmax = max(a)

# step = (valmax-valmin)/40


# while valmin<valmax:
#     ab1mask = abdf['val']>valmin
#     ab2mask = abdf['val']<valmin+step

#     c.append(np.mean([valmin, valmin+step]))
#     abdf_k = abdf[ab1mask&ab2mask]

#     d1.append(abdf_k['err'].quantile(0.05))
#     d.append(abdf_k['err'].mean())
#     valmin+=step

# plt.plot(c,d)
# plt.xlabel('LFQ log2')
# plt.ylabel('avg %missing vals')
# plt.title('avg %missing vals by lfq')
# plt.savefig(path+'depict\\protMISS_02.png',dpi=300)
# plt.show()
# plt.clf()

# plt.plot(c,d1)
# plt.xlabel('LFQ log2')
# plt.ylabel('95 %missing vals')
# plt.title('95 %missing vals by lfq')
# plt.savefig(path+'depict\\protMISS_03.png',dpi=300)
# plt.show()
# plt.clf()

# plt.plot(c,[200/x for x in d1])
# plt.xlabel('LFQ log2')
# plt.ylabel('#reps to 2xID w 95%conf')
# plt.title('#reps to 2xID w 95%conf lfq')
# plt.savefig(path+'depict\\protMISS_04.png',dpi=300)
# plt.show()
# plt.clf()


# plt.plot(c,[100/x for x in d1])
# plt.xlabel('LFQ log2')
# plt.ylabel('#reps to 1xID w 95%conf')
# plt.title('#reps to 1xID w 95%conf lfq')
# plt.savefig(path+'depict\\protMISS_05.png',dpi=300)
# plt.show()
# plt.clf()

# pd.DataFrame({'bin':c, 'avID':d}).to_csv(path+'hitrate.csv',index=False)



# PCA analysis =======================================================================================================


# protval_dn = protval.dropna().reset_index(drop=True)
# protval_dn['mean'] = protval_dn.mean(axis=1)
# protval_dn = protval_dn[protval_dn['mean']>2**14].reset_index(drop=True)
# del protval_dn['mean']

# pvdn = np.log2(protval_dn).values


# names = ['1', '2', '3', '4', '5', '6',
#           '7', '8', '9', '10', '11', '12',
#           '13', '14', '15', '16', '17', '18', '19', '20']
# ncolors = []

# colorchoices = [colorsplit(col/(len(protval.columns)-1)) for col in range(len(protval.columns))]

# color_count = 0
# color_dict = {}

# for name in names:
#     if ptna_name_rep(name) in color_dict.keys():
#         ncolors.append(color_dict[ptna_name_rep(name)])
#     else:
#         color_dict[ptna_name_rep(name)] = colorchoices[color_count]
#         ncolors.append(colorchoices[color_count])
#         color_count += 1

# pca = PCA(n_components=2)
# pCones = pca.fit_transform(pvdn.T)
# principalDf = pd.DataFrame(data=pCones, columns = ['PC1', 'PC2'])

# xs = principalDf['PC1']
# ys = principalDf['PC2']

# plt.scatter(xs, ys, c=np.array(ncolors))

# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')

# plt.legend(handles=[mpatches.Patch(color=col, label=text) for text, col in color_dict.items()],
#                     prop={'size': 6})

# plt.savefig(path+'depict\\PCAmap.png',dpi=300)
# plt.show()
# plt.clf()



# a = []
# b = []

# for x in range(len(xs)):
#     a.append(x+1)
#     b.append(xs[x]**2 + ys[x]**2)

# plt.scatter(a, b)
# plt.xlabel('rep #')
# plt.ylabel('PCA condensed CoV')
# plt.title('Condensed CoV over time')
# plt.savefig(path+'depict\\pcaTIME.png',dpi=300)
# plt.show()
# plt.clf()

# pca = PCA(n_components=1)
# pCones = pca.fit_transform(pvdn.T)
# principalDf = pd.DataFrame(data=pCones, columns = ['PC1'])

# xs = principalDf['PC1']


# plt.scatter(a, xs)
# plt.xlabel('rep #')
# plt.ylabel('1D PCA')
# plt.title('1D PCA over time')
# plt.savefig(path+'depict\\pca1D.png',dpi=300)
# plt.show()
# plt.clf()




# protval_dn['mean'] = protval_dn.mean(axis=1)





# IMPUTE strategy 3x rep  ===================================================================================================



# mask = protval.isna().any(axis=1)
# pvdn = protval.loc[mask].reset_index(drop=True)

# a=[]
# # b=[]

# for i in pvdn.index:
#     protslice = pvdn.iloc[i]
#     mstat = protslice.dropna().mean()
#     if len(protslice.dropna())>1:
#         for _ in range(len(protslice)*2):
#             rs = pd.Series(random.sample([x for x in protslice.values], 3)).dropna()
#             if len(rs)>0:
#                 ustat = np.mean(rs)
#                 a.append(100*(abs(mstat-ustat)/ustat))

# plt.hist(a, alpha=0.5, color='blue', bins=30, range=(0,100))
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: dropna')
# plt.savefig(path+'depict\\err_dropna.png',dpi=300)
# plt.show()
# plt.clf()


# mask = protval.isna().any(axis=1)
# pvdn = protval.loc[mask].reset_index(drop=True)

# # a=[]
# b=[]

# for i in pvdn.index:
#     protslice = pvdn.iloc[i]
#     mstat = protslice.dropna().mean()
#     if len(protslice.dropna())>1:
#         for _ in range(len(protslice)*2):
#             rs = pd.Series(random.sample([x for x in protslice.values], 3)).fillna(0)
#             ustat = np.mean(rs.dropna())
#             b.append(100*(abs(mstat-ustat)/ustat))

# plt.hist(b, alpha=0.5, color='red', bins=30, range=(0,100))
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: fillna=0')
# plt.savefig(path+'depict\\err_fillna0.png',dpi=300)
# plt.show()
# plt.clf()



# # mask = protval.isna().any(axis=1)
# # pvdn0 = protval.loc[mask]


# # slot= 10
# # protean=protval.copy()

# # while protean.isna().any().any():
# #     # print(protean.isna().any())
# #     # print('-',end='')
# #     imp_col = random.sample([x for x in protean.columns], 3)
# #     impdf = knn_tn_impute(protean[imp_col], slot)
# #     for x,col in zip(impdf.T, imp_col):
# #         protean[col] = x

# # pvdn = protean.loc[pvdn0.index].reset_index()

# # c=[]

# # for i in pvdn.index:
# #     protslice = pvdn.iloc[i]
# #     mstat = protslice.dropna().mean()
# #     if len(protslice.dropna())>1:
# #         for _ in range(len(protslice)*2):
# #             rs = pd.Series(random.sample([x for x in protslice.values], 3)).dropna()
# #             ustat = np.mean(rs)
# #             c.append(100*(abs(mstat-ustat)/ustat))
# # print(f"{slot}: {np.median(c)}")

# # plt.hist(c, alpha=0.5, color='green', bins=30, range=(0,100))
# # plt.ylabel('#peptide permutations')
# # plt.xlabel('%error')
# # plt.title('Error Distribution: knntn ')
# # plt.savefig(path+f'depict\\err_KNNTNna{slot}.png',dpi=300)
# # plt.show()
# # plt.clf()


# mask = protval.isna().any(axis=1)
# pvdn0 = protval.loc[mask]



# slot= 10
# protean=protval.copy()

# while protean.isna().any().any():
#     # print(protean.isna().any())
#     print('-',end='')
#     imp_col = random.sample([x for x in protean.columns], len(protean.columns))
#     impdf = knn_tn_impute(protean[imp_col], slot)
#     for x,col in zip(impdf.T, imp_col):
#         protean[col] = x

# pvdn = protean.loc[pvdn0.index].reset_index()

# d=[]
# for i in pvdn.index:
#     protslice = pvdn.iloc[i]
#     mstat = protslice.dropna().mean()
#     if len(protslice.dropna())>1:
#         for _ in range(len(protslice)*2):
#             rs = pd.Series(random.sample([x for x in protslice.values], 3)).dropna()
#             ustat = np.mean(rs)
#             d.append(100*(abs(mstat-ustat)/ustat))
# print(f"{np.median(d)}")

# plt.hist(d, alpha=0.5, color='green', bins=30, range=(0,100))
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: knntn')
# plt.savefig(path+f'depict\\err_KNNTNna{slot}.png',dpi=300)
# plt.show()
# plt.clf()


# ascore = np.median(a)
# bscore = np.median(b)
# cscore = np.median(d)


# plt.hist(a, color='blue',alpha=0.5, label='dropna', bins=30, range=(0,100))
# plt.hist(b, color='red', alpha=0.5, label='fillna0', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=ascore, color='blue')
# plt.axvline(x=bscore, color='red')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: fill vs drop')
# plt.savefig(path+f'depict\\err_comp_01.png',dpi=300)
# plt.show()
# plt.clf()


# plt.hist(a, color='blue',alpha=0.5, label='dropna', bins=30, range=(0,100))
# plt.hist(c, color='green', alpha=0.5, label='5NN-TN', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=ascore, color='blue')
# plt.axvline(x=cscore, color='green')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: 5NN vs drop')
# plt.savefig(path+f'depict\\err_comp_02.png',dpi=300)
# plt.show()
# plt.clf()


# plt.hist(b, color='red', alpha=0.5, label='fillna0', bins=30, range=(0,100))
# plt.hist(c, color='green', alpha=0.5, label='5NN-TN', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=cscore, color='green')
# plt.axvline(x=bscore, color='red')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: 5NN vs drop')
# plt.savefig(path+f'depict\\err_comp_03.png',dpi=300)
# plt.show()
# plt.clf()




# IMPUTE strategy again  ===================================================================================================


# pvdn = protval.copy()
# a,b,c = [],[],[]

# pca = PCA(n_components=6)
# pCones = pca.fit_transform(pvdn.dropna().T)
# principalDf = pd.DataFrame(data=pCones, columns = ['PC1', 'PC2', 'PC3','PC4', 'PC5', 'PC6'], index=pvdn.columns)
# exvar = pca.explained_variance_ratio_

# k = 0
# for i in pvdn.index:
#     if int(100*i/pvdn.shape[0])>int(100*k/pvdn.shape[0]):
#         print(f"{int(100*i/pvdn.shape[0])}%")
#     protslice = pvdn.iloc[i]
#     mstat = protslice.dropna().mean()
#     if len(protslice.dropna())>1:
#         for _ in range(len(protslice.dropna())*2):
#             ls = random.sample([x for x in protslice.index], 3)
#             rs = pd.Series(protslice[ls])
#             if len(rs.dropna())==2:
#                 a_ustat = np.mean(rs.dropna())
#                 b_ustat = np.mean(rs.fillna(0))
#                 a.append(100*(abs(mstat-a_ustat)/a_ustat))
#                 b.append(100*(abs(mstat-b_ustat)/b_ustat))
#                 rs = phuff_impute(rs, ls, principalDf.loc[ls], exvar)
#                 c_ustat = np.mean(rs)
#                 c.append(100*(abs(mstat-c_ustat)/c_ustat))
#     k = i


# ascore = np.median(a)
# bscore = np.median(b)
# cscore = np.median(c)

# plt.hist(a, color='blue',alpha=0.5, label='dropna', bins=30, range=(0,100))
# plt.hist(b, color='red', alpha=0.5, label='fillna0', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=ascore, color='blue')
# plt.axvline(x=bscore, color='red')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: fill vs drop')
# plt.savefig(path+f'depict\\err_phuff_01.png',dpi=300)
# plt.show()
# plt.clf()


# plt.hist(a, color='blue',alpha=0.5, label='dropna', bins=30, range=(0,100))
# plt.hist(c, color='green', alpha=0.5, label='3xREP iso', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=ascore, color='blue')
# plt.axvline(x=cscore, color='green')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: 3xREP iso vs drop')
# plt.savefig(path+f'depict\\err_phuff_02.png',dpi=300)
# plt.show()
# plt.clf()


# plt.hist(b, color='red', alpha=0.5, label='fillna0', bins=30, range=(0,100))
# plt.hist(c, color='green', alpha=0.5, label='3xREP iso', bins=30, range=(0,100))
# plt.legend()
# plt.axvline(x=cscore, color='green')
# plt.axvline(x=bscore, color='red')
# plt.ylabel('#peptide permutations')
# plt.xlabel('%error')
# plt.title('Error Distribution: 3xREP iso vs fill')
# plt.savefig(path+f'depict\\err_phuff_03.png',dpi=300)
# plt.show()
# plt.clf()



# INSM1 targeting           ============================================================================

# protean = prot.copy()

# poidf = pd.read_csv(path+'pois_illus.csv')

# targets = poidf['accession'].values

# protean['target'] = [sum([x.count(y) for y in targets]) for x in protean['Protein.Ids']]

# print(protean.shape)
# protean = protean[protean['target']>0]
# print(protean.shape)

# protean.to_csv(path+'pois_targeted1.csv',index=False)



# r2 of proteome mix 1:1 =============================================================================================

c = []
for _ in range(3*len(protval.columns)):
    cols = random.sample([x for x in protval.columns], 2)
    pvdn = protval[cols]

    a = []
    b = []
    for i in pvdn.index:
        protslice = pvdn.iloc[i].dropna()
        if len(protslice)==2:
            rs = [x for x in protslice.values]
            a.append(np.log2(rs[0]))
            b.append(np.log2(rs[1]))

    x = np.asarray(a).reshape(-1,1)
    y = np.asarray(b).reshape(-1,1)

    regr = linear_model.LinearRegression()
    regr.fit(x, y)
    score = r2_score(y, regr.predict(x))
    print(score)
    c.append(score)
print()
print(np.median(c))

plt.hist(c)
plt.suptitle('linear r2 corr of 1:1 replicates')
plt.title('60 random comparisons total')
plt.ylabel('# comparisons')
plt.xlabel('r2 score')
plt.savefig('corrscores_dia.png')
plt.show()
plt.clf()






















