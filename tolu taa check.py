# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:49:18 2024

@author: HUFFMP
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
from sklearn.svm import SVR
from scipy.stats import ttest_ind as tti










df = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/procan TAA vec_proCanno.xlsx')
tf = pd.read_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu TAA true.xlsx')

df = df.fillna(0)

CLx = ['Calu-6', 'ChaGo-K-1', 'HCC-78', 'NCI-H1648', 'NCI-H1650', 'SK-MES-1']

nf = tf[tf.columns[0]]
tf = tf[tf.columns[1:]]
# tf=np.log2(tf)
# tf = tf.replace(-np.inf, np.nan)
# tf = tf.fillna(tf.min().min()-1)

r_r = []
p_r = []
model = LinearRegression()

truX = []
truY1 = []
truY2 = []
truCL = []
truPROT = []

for CL in CLx:
    X = tf[CL].values.reshape(-1, 1)
    y1 = df[f'ProCAN Abundance: {CL}'].values
    y2 = df[f"ProCAN TPM: {CL}"].values
    truPROT += [x for x in nf.values]

    truX += [x for x in X]
    truY1 += [y for y in y1]
    truY2 += [y for y in y2]
    truCL += [CL for y in y2]

    model = LinearRegression()
    # model = SVR()
    model.fit(X, y1)
    r2_1 = r2_score(y1, model.predict(X))

    model = LinearRegression()
    # model = SVR()
    model.fit(X, y2)
    r2_2 = r2_score(y2, model.predict(X))

    r_r.append(r2_2)
    p_r.append(r2_1)

    # plt.scatter(X, y1)
    # plt.title(f'{CL}: RQuant v PTMX Abundance   -   r2:{round(r2_1,2)}')
    # plt.ylabel('Abundance')
    # plt.xlabel('Rquant lin')
    # # plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu taa/1-{CL}_ptmx.png')
    # plt.show()
    # plt.scatter(X, y2)
    # plt.title(f'{CL}: RQuant v TPM   -   r2:{round(r2_2,2)}')
    # plt.ylabel('TPM')
    # plt.xlabel('Rquant lin')
    # # plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu taa/1-{CL}_tpm.png')
    # plt.show()

# print('r: ', np.mean(r_r))
# print('p: ', np.mean(p_r))


model = LinearRegression()
# model = SVR()
model.fit(truX, truY1)
r2_1 = r2_score(truY1, model.predict(truX))

model = LinearRegression()
# model = SVR()
model.fit(truX, truY2)
r2_2 = r2_score(truY2, model.predict(truX))

# plt.scatter(truX, truY1)
# plt.title(f'ALL: RQuant v PTMX Abundance   -   r2:{round(r2_1,2)}')
# plt.ylabel('Abundance')
# plt.xlabel('Rquant lin')
# plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu taa/1-ALL_ptmx.png')
# plt.show()
plt.scatter(np.log10(truY2), np.log10(truX))
plt.title(f'ALL: RQuant v TPM   -   r2:{round(r2_2,2)}')
plt.ylabel('TPM')
plt.xlabel('Rquant lin')
plt.savefig(f'//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tolu taa/1-ALL_tpm.png')
plt.show()

df = pd.DataFrame({'RQuant':[x[0] for x in np.log10(truX)], 'PTMX':truY1, 'TPM':np.log10(truY2), 'Cell':truCL, 'Prot':truPROT})
alpha = 1000

df = df[df['RQuant']!=-np.inf]
df = df[df['TPM']!=-np.inf]

for c in df['Cell'].unique():
    x = df[df['Cell']==c]['RQuant'].values.reshape(-1,1)
    y = df[df['Cell']==c]['TPM'].values
    model = LinearRegression()
    # model = SVR()
    model.fit(x, y)
    r2 = r2_score(y, model.predict(x))
    print(f'{c}: {r2}')
    plt.scatter(y,x,label=c)
plt.legend()
plt.ylabel('log10 TPM')
plt.xlabel('log10 RQuant')
plt.title(f'cell cluster')
plt.savefig('TAA_merged_cells.png', dpi=400)

plt.show()


# NUM_COLORS = len(df['Prot'].unique())+3
# cm = plt.get_cmap('gist_rainbow')
# color = [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)]
# color1 = [[i/1.3 for i in c][:-1]+[1] for c in color]


ci = 0
for c in df['Prot'].unique():
    x = df[df['Prot']==c]['RQuant'].values.reshape(-1,1)
    y = df[df['Prot']==c]['TPM'].values
    model = LinearRegression()
    # model = SVR()
    model.fit(x, y)
    r2 = r2_score(y, model.predict(x))
    print(f'{c}: {r2}')
    plt.scatter(y,x,label=c)
    ci+=1
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.ylabel('log10 TPM')
plt.xlabel('log10 RQuant')
plt.title('prot cluster')
plt.savefig('TAA_merged_proteins.png', dpi=400, bbox_inches='tight')

plt.show()


# quanta = 0.20

# print(f'Confidence Quantile: {quanta}')
# print('ptmx cutoffs')
# g = 0
# i = 0
# ddf = df.copy()
# while g<1000:
#     g = ddf[ddf['PTMX']>(i*df['PTMX'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 1000: PTMX cutoff of {i*df["PTMX"].max()/alpha}')

# g = 0
# i = 0
# ddf = df.copy()
# while g<5000:
#     g = ddf[ddf['PTMX']>(i*df['PTMX'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 5000: PTMX cutoff of {i*df["PTMX"].max()/alpha}')


# g = 0
# i = 0
# ddf = df.copy()
# while g<10000:
#     g = ddf[ddf['PTMX']>(i*df['PTMX'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 10000: PTMX cutoff of {i*df["PTMX"].max()/alpha}')



# print()

# print('tpm cutoffs')
# g = 0
# i = 0
# ddf = df.copy()
# while g<1000:
#     g = ddf[ddf['TPM']>(i*df['TPM'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 1000: TPM cutoff of {i*df["TPM"].max()/alpha}')

# g = 0
# i = 0
# ddf = df.copy()
# while g<5000:
#     g = ddf[ddf['TPM']>(i*df['TPM'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 5000: TPM cutoff of {i*df["TPM"].max()/alpha}')


# g = 0
# i = 0
# ddf = df.copy()
# while g<10000:
#     g = ddf[ddf['TPM']>(i*df['TPM'].max()/alpha)]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 10000: TPM cutoff of {i*df["TPM"].max()/alpha}')

# print()



# m3 = 2
# print('dual cutoffs')
# g = 0
# i = 0
# ddf = df.copy()
# while g<1000:
#     m1 = ddf['TPM']>(i*df['TPM'].max()/alpha)
#     m2 = ddf['PTMX']>(m3*i*df['PTMX'].max()/alpha)
#     g = ddf[m1&m2]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 1000: TPM cutoff of {i*df["TPM"].max()/alpha}, PTMX cutoff of {m3*i*df["PTMX"].max()/alpha}')

# g = 0
# i = 0
# ddf = df.copy()
# while g<5000:
#     m1 = ddf['TPM']>(i*df['TPM'].max()/alpha)
#     m2 = ddf['PTMX']>(m3*i*df['PTMX'].max()/alpha)
#     g = ddf[m1&m2]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 5000: TPM cutoff of {i*df["TPM"].max()/alpha}, PTMX cutoff of {m3*i*df["PTMX"].max()/alpha}')


# g = 0
# i = 0
# ddf = df.copy()
# while g<10000:
#     m1 = ddf['TPM']>(i*df['TPM'].max()/alpha)
#     m2 = ddf['PTMX']>(m3*i*df['PTMX'].max()/alpha)
#     g = ddf[m1&m2]['RQuant'].quantile(quanta)
#     i+= 1
# i-=1
# print(f'RQ 10000: TPM cutoff of {i*df["TPM"].max()/alpha}, PTMX cutoff of {m3*i*df["PTMX"].max()/alpha}')


# FDR ==========================================================================================

# ptdf = pd.DataFrame()
# for ig in [(x+1)/40 for x in range(19)]:
#     quanta = ig
#     stat = []
#     for j in [(x+1)*500 for x in range(20)]:
#         g = 0
#         i = 0
#         ddf = df.copy()
#         while g<j:
#             g = ddf[ddf['PTMX']>(i*df['PTMX'].max()/alpha)]['RQuant'].quantile(quanta)
#             i+= 1
#         i-=1
#         stat.append(i*df['PTMX'].max()/alpha)
#     ptdf[quanta] = stat
# ptdf.index = [(x+1)*500 for x in range(20)]
# # ptdf.to_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/lfq atom.xlsx')



# tpdf = pd.DataFrame()
# for ig in [(x+1)/40 for x in range(19)]:
#     quanta = ig
#     stat = []
#     for j in [(x+1)*500 for x in range(20)]:
#         g = 0
#         i = 0
#         ddf = df.copy()
#         while g<j:
#             g = ddf[ddf['TPM']>(i*df['TPM'].max()/alpha)]['RQuant'].quantile(quanta)
#             i+= 1
#         i-=1
#         stat.append(i*df["TPM"].max()/alpha)
#     tpdf[quanta] = stat
# tpdf.index = ptdf.index = [(x+1)*500 for x in range(20)]
# # tpdf.to_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/tpm atom.xlsx')


# ==============================================

# df['mstat'] = df['RQuant'] / (df['PTMX']+1)

# protl = df['Prot'].unique()

# p1 = []
# p2 = []
# logp = []
# diff = []
# p1d = []
# p2d = []

# for p in protl:
#     m1 = df[df['Prot']==p]['mstat'].values
#     for pi in protl:
#         m2 = df[df['Prot']==pi]['mstat'].values
#         p1.append(p)
#         p2.append(pi)
#         logp.append(-np.log10(tti(m1, m2)[1]))
#         diff.append((1+np.mean(m2))/(1+np.mean(m1)))
#         p1d.append(np.mean(m1))
#         p2d.append(np.mean(m2))

# fd = pd.DataFrame({'prot1':p1, 'val1':p1d, 'prot2':p2, 'val2':p2d, 'nlogp':logp, 'diff':diff})
# fd.to_excel('//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/HUFFMP/misc/mstat comp.xlsx')






















