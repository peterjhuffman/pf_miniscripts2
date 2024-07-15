# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:11:06 2024

@author: HUFFMP
"""
import pandas as pd
from datetime import date


def parenth_cut(stri):
    b = ''
    sw = True
    for char in stri:
        if sw:
            b+=char
        if char=='(':
            sw=False
        if char==')':
            sw=True
            b = b[:-2]
    return b

def quant_interp(val):
    if val =='NA%':
        return '0'
    if val =='N/A':
        return '0'
    else:
        return val[:-1]

filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\meth\\'
filename = 'wjxk_2402.xlsx'

todate = str(date.today()).split('-')

n = input('TMT16plex code name: ')
v = input('TMT16plex input string: ')

n2 = input('TMT18plex code name: ')
v2 = input('TMT18plex input string: ')

tmtname = f"TMTpro18plex_{n}_{n2}_{todate[0][2:]}{todate[1]}.method"

# v = 'TMTpro126 126.127726 N/A N/A N/A N/A 100% 0.31% 9.09% (127C) 0.02% 0.32% TMTpro127N 127.124761 N/A N/A N/A 0.57% (126) 100% N/A 9.79% (128N) N/A 0.33% TMTpro127C 127.131081 N/A N/A 0.84% (126) N/A 100% 0.23% 8.40% (128C) 0.02% 0.27% TMTpro128N 128.128116 N/A 0.00% 0.68% (127N) 0.52% 100% N/A 5.23% (129N) N/A 0.20% TMTpro128C 128.134436 0.00% N/A 1.44% (127C) N/A 100% 0.34% 6.26% (129C) 0.00% 0.17% TMTpro129N 129.131471 0.00% 0.14% 1.30% (128N) 0.89% 100% N/A 7.52% (130N) N/A 0.12% TMTpro129C 129.13779 0.13% N/A 2.59% (128C) N/A 100% 0.32% 6.07% (130C) 0.01% 0.09% TMTpro130N 130.134825 0.13% 0.00% 2.41% (129N) 0.27% 100% N/A 5.58% (131N) N/A 0.10% TMTpro130C 130.141145 0.25% N/A 3.22% (129C) N/A 100% 0.28% 5.06% (131C) 0.00% 0.06% TMTpro131N 131.13818 0.04% 0.01% 2.73% (130N) 0.49% 100% N/A 3.13% (132N) N/A 0.06% TMTpro131C 131.1445 0.09% N/A 4.02% (130C) N/A 100% 1.17% 3.62% (132C) 0.02% 0.03% TMTpro132N 132.141535 0.07% 0.01% 3.14% (131N) 0.73% 100% N/A 3.40% (133N) N/A 0.03% TMTpro132C 132.147855 0.13% N/A 5.14% (131C) N/A 100% 1.16% 1.92% (133C) 0.00% 0.00% TMTpro133N 133.14489 0.14% 0.02% 3.37% (132N) 0.63% 100% N/A 1.16% (134N) N/A 0.00% TMTpro133C 133.15121 0.18% N/A 4.14% (132C) N/A 100% 0.40% 1.11% (134C) 0.00% N/A TMTpro134N 134.148245 0.28% 0.10% 5.52% (133N) 0.35% 100% N/A 1.12% (135N) N/A N/A'
# v2 = 'TMTpro-134C 134.154565 0.16% N/A 5.81% (133C) N/A 100% 0.39% N/A (135C) N/A N/A TMTpro-135N 135.151600 0.21% 0.04% 5.90% (134N) 0.74% 100% N/A N/A (136N) N/A N/A'

# x = input('TMT String:')
x=v
cstring = x.split('TMTpro1')[1:]
cf = pd.DataFrame([parenth_cut(x).split(' ')[:13] for x in cstring])

cf.index = list(cf[0])
cf.columns = ['tag', 'mass', '-2x 13c', '-13c-15n', '-13c','-15n','-','+15n','+13c','+15n+13c','+2x13c', '']
del cf['']
del cf['tag']

# y = input('TMT String (17,18):')
y=v2
tstring = y.split('TMTpro-1')[1:]
tf = pd.DataFrame([parenth_cut(x).split(' ')[:13] for x in tstring])
tf.index = list(tf[0])
tf.columns = ['tag', 'mass', '-2x 13c', '-13c-15n', '-13c','-15n','-','+15n','+13c','+15n+13c','+2x13c','']
del tf['']
del tf['tag']

# cf.to_excel(filepath+filename, index=False)


file1 = open(filepath+'TMTpro 18plex.txt', 'r')
Lines = file1.readlines()
file1.close()

# print([x for x in Lines])
newLines = []

for val, line in enumerate(Lines):
    newLines.append(f"{val}: {line}")


tag = '26'
cols = ['-2x 13c', '-13c-15n', '-13c','-15n','+15n','+13c','+15n+13c','+2x13c']
xs = 14


Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '27N'
xs = 76
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '27C'
xs = 138
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '28N'
xs = 200
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '28C'
xs = 262
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '29N'
xs = 324
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '29C'
xs = 386
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]


tag = '30N'
xs = 448
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '30C'
xs = 510
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '31N'
xs = 572
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '31C'
xs = 634
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '32N'
xs = 696
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '32C'
xs = 758
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '33N'
xs = 820
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '33C'
xs = 882
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '34N'
xs = 944
Lines[xs] = Lines[xs][:35]+quant_interp(cf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(cf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(cf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(cf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(cf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(cf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(cf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(cf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '34C'
xs = 1006
Lines[xs] = Lines[xs][:35]+quant_interp(tf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(tf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(tf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(tf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(tf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(tf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(tf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(tf.loc[tag,cols[7]])+Lines[xs+48][36:]

tag = '35N'
xs = 1068
Lines[xs] = Lines[xs][:35]+quant_interp(tf.loc[tag,cols[0]])+Lines[xs][36:]
Lines[xs+6] = Lines[xs+6][:35]+quant_interp(tf.loc[tag,cols[1]])+Lines[xs+6][36:]
Lines[xs+12] = Lines[xs+12][:35]+quant_interp(tf.loc[tag,cols[2]])+Lines[xs+12][36:]
Lines[xs+18] = Lines[xs+18][:35]+quant_interp(tf.loc[tag,cols[3]])+Lines[xs+18][36:]
Lines[xs+30] = Lines[xs+30][:35]+quant_interp(tf.loc[tag,cols[4]])+Lines[xs+30][36:]
Lines[xs+36] = Lines[xs+36][:35]+quant_interp(tf.loc[tag,cols[5]])+Lines[xs+36][36:]
Lines[xs+42] = Lines[xs+42][:35]+quant_interp(tf.loc[tag,cols[6]])+Lines[xs+42][36:]
Lines[xs+48] = Lines[xs+48][:35]+quant_interp(tf.loc[tag,cols[7]])+Lines[xs+48][36:]


writepath = '//amer.pfizer.com/pfizerfiles/Research/LAJ/oru-tcb/ChemBio/general - chemoproteomics/TMT reagentBatches/'
file2 = open(writepath+tmtname, 'a')
for nL in Lines:
    file2.write(nL)
file2.close()