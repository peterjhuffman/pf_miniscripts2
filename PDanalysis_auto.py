print('Making graphs. Analyzing. Whatnot.')

import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from sklearn.decomposition import PCA
import matplotlib.patches as mpatches
pd.options.mode.chained_assignment = None


def proteonalysis_main(filename, df):
    df = df.dropna(subset=['126'])

    xdf = ptna_normalize(df)

    names = ptna_load_names(filename)
    if names == False:
        names = ptna_naming()
        ptna_save_names(filename, names)
#    ptna_pca_eval(xdf, names, filename[:-4])
    df = ptna_name_assignment(df, names)
#    ptna_vis_target(df, names, filename[:-4])
    ptna_vis_heats(df, names, filename[:-4])


def ptna_load_names(filename):
    file_present = os.path.isfile(f'PDanalysis_results{filename[:-4]}/tmt_lane_names.txt')

    if file_present:
        direct = print('There is a name file loaded.')
        direct = 'y'
        while (direct.upper() != 'N') and (direct.upper() != 'Y'):
            print('Unknown response. Try again.')
            direct = input('There is a name file loaded. Would you like to use? (Y/N):')
        if direct.upper() == 'Y':
            f = open(f'PDanalysis_results{filename[:-4]}/tmt_lane_names.txt', 'r')
            details = f.read().split(', ')
            f.close()
            print(f'Names loaded as: {details}')
            return details
    return False


def ptna_normalize(df, lanes=0):
    if type(lanes) == int:
        lanes = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C', '131N', '131C',
                 '132N', '132C', '133N', '133C', '134C', '134N', '135N']
    xdf = df[lanes]

    xdf = (xdf - xdf.median(axis=0)) / xdf.std(axis=0)

#    for lane in lanes:
#        df[lane] = xdf[lane]
    return xdf


def ptna_naming():
    names = []
    for x in range(18):
        names.append(input(f'TMT Lane {x+1}: '))

    names = ptna_name_check(names)


    namecounts = {}
    for x in names:
        namecounts[x] = 0

    for x in range(len(names)):
        namecounts[names[x]] += 1
        names[x] = f"{names[x]}_{namecounts[names[x]]}"#[:-2]

    return names


def ptna_name_check(names):
    namecheck = input('Are those the correct names? (Y/N): ')
    if (namecheck.upper() != 'N') and (namecheck.upper() != 'Y'):
        print('Unknown response. Try again.')
        names = ptna_name_check(names)
    elif namecheck.upper() == 'N':
        print('Enter the names again.')
        return ptna_naming()
    return names


def ptna_vis_heats(df, names, folder):
    col_df = ptna_collate_conditions(df)
    log_col_df = pd.DataFrame(columns=col_df.columns)
    print(col_df.head())

    col_df.to_csv('bread.csv')

    log_col_df['24'] = np.log2(col_df['24']/col_df['24'])
    log_col_df['48'] = np.log2(col_df['48']/col_df['24'])
    log_col_df['5day'] = np.log2(col_df['5day']/col_df['24'])
    log_col_df['24dox'] = np.log2(col_df['24dox']/col_df['24dox'])
    log_col_df['48dox'] = np.log2(col_df['48dox']/col_df['24dox'])
    log_col_df['5daydox'] = np.log2(col_df['5daydox']/col_df['24dox'])
    log_col_df.index=df['Accession']
#    ptna_normalize(log_col_df, log_col_df.columns).to_csv('bread.csv')
    log_col_df.to_csv('bread2.csv')

#    log_col_df = ptna_normalize(np.log2(col_df), col_df.columns)
#    plt.imshow(log_col_df.T, cmap='magma', interpolation='nearest', aspect=1000)
#    plt.xticks([])
#    locs, labels = plt.yticks()
#    plt.yticks(np.arange(len(log_col_df.columns)), list(log_col_df.columns))
#    plt.savefig(f'PDanalysis_results\\{folder}\\'+'proteomewide_comp.png', dpi=2000)
#
#    varia = ['dox']
#    for var in varia:
#        col_slice = col_df[var].values
#        control_slice = col_df['control'].values


def ptna_save_names(filename, names):
    print(f"Names: {names}")
    save_names_q = input('Would you like to save these lane names for future runs? (Y/N): ')
    while (save_names_q.upper() != 'Y') and (save_names_q.upper() != 'N'):
        print('Unknown response. Try again.')
        save_names_q = input('Would you like to save these lane names for future runs? (Y/N): ')
    if save_names_q.upper() == 'Y':
        with open(f'PDanalysis_results{filename[:-4]}/tmt_lane_names.txt', 'w') as f:
            f.write(', '.join(names))


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


def ptna_pca_eval(df, names, folder):
    df = df.dropna(subset=['126'])
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

    start = df.columns.get_loc('126')
    end = df.columns.get_loc('135N')+1

    x = df.iloc[:, start:end].values
    df.iloc[:, start:end].to_csv('pca_input.csv', index=False)
    pca = PCA(n_components=2)

    pCones = pca.fit_transform(x.T)
    principalDf = pd.DataFrame(data=pCones, columns = ['PC1', 'PC2'])
    principalDf.to_csv('principalCA.csv', index=False)

    xs = principalDf['PC1']
    ys = principalDf['PC2']

    ncolors = pd.Series(ncolors)
    for color in ncolors.unique():
        cmap = ncolors==color
        edges = ((pCones[cmap][:, 0].min()),
                 (pCones[cmap][:, 0].max()),
                 (pCones[cmap][:, 1].min()),
                 (pCones[cmap][:, 1].max()))
        if (edges[1]-edges[0]) > (edges[3]-edges[2]):
            linearity = 0
        else:
            linearity = 1

        rempop = [0, 1, 2, 3, 4, 5, 6, 7, 8]

        pointmin = pCones[cmap][:, linearity].argmin()
        pointmax = pCones[cmap][:, linearity].argmax()
        rempop.remove(pointmin)
        rempop.remove(pointmax)
        pointrem = rempop[0]

        pointmin = pCones[cmap][pointmin]
        pointmax = pCones[cmap][pointmax]
        pointdist = math.sqrt((pointmax[0]-pointmin[0])**2+(pointmax[1]-pointmin[1])**2)

        center = ((pointmax[0]+pointmin[0])/2, (pointmax[1]+pointmin[1])/2)
        angle = math.atan((pointmax[0]-pointmin[0]) / (pointmax[1]-pointmin[1])) * (180/math.pi)

        tridist = ptna_ellipsewidth(pCones[cmap][pointrem], (pointmin, pointmax))
        circmetric = tridist/pointdist


        if circmetric < 0.333333:
            circle = mpl.patches.Ellipse(center, (pointdist)+1, tridist+1, 90-angle, fc=(0, 0, 0, 0),ec=color)
        else:
            center = ((pCones[cmap][:, 0].sum() / 3), (pCones[cmap][:, 1].sum() / 3))
            size = 0
            for pt in pCones[cmap]:
                distc = math.sqrt((pt[0]-center[0])**2+(pt[1]-center[1])**2)
                size = max(size, distc)
            circle = mpl.patches.Circle(center, size+0.7, fc=(0,0,0,0), ec=color)
        plt.gca().add_patch(circle)

    plt.scatter(xs, ys, c=np.array(ncolors))

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA of Lane Similarity')

    plt.legend(handles=[mpatches.Patch(color=col, label=text) for text, col in color_dict.items()],
                        prop={'size': 6})
    plt.savefig(f'PDanalysis_results\\{folder}'+'\\TMT_pcagroups.png', dpi=1200)
    plt.clf()
    plt.close()

    plt.scatter(xs, ys, c=np.array(ncolors))

    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')

    plt.legend(handles=[mpatches.Patch(color=col, label=text) for text, col in color_dict.items()],
                        prop={'size': 6})
    plt.title('PCA of Lane Similarity')
    plt.savefig(f'PDanalysis_results\\{folder}'+'\\TMT_pca.png', dpi=1200)
    plt.clf()
    plt.close()


def ptna_name_assignment(df, names):
    cols = pd.Series(df.columns)
    locus = cols[cols=='126'].index[0]

    for x in range(18):
        df.rename(columns = {df.columns[locus+x] : names[x]}, inplace=True)

    return df


def ptna_underclip(name):
    while name.count('_') > 0:
        name = name[:-1]
    return name


def ptna_collate_conditions(df, coldict=0):
    if coldict == 0:
#        coldict = {'control':['24_1', '24_2', '24_3', '48_1', '48_2', '5day_1',
#                              '48_3', '5day_2', '5day_3'],
#                   'dox':['24dox_1', '24dox_2', '48dox_1', '24dox_3', '48dox_2',
#                          '48dox_3', '5daydox_1', '5daydox_2', '5daydox_3']}
        coldict = {'24':['24_1', '24_2', '24_3'],
                   '48': ['48_1', '48_2', '48_3'],
                   '5day':['5day_1', '5day_2', '5day_3'],
                   '24dox':['24dox_1', '24dox_2', '24dox_3'],
                   '48dox':['48dox_1', '48dox_2', '48dox_3'],
                   '5daydox':['5daydox_1', '5daydox_2', '5daydox_3']}
    newdf = pd.DataFrame(columns=list(coldict.keys()))
    for key in coldict.keys():
        newdf[key] = df[coldict[key]].sum(axis=1)
    return newdf


def ptna_vis_targets(df, names, folder):
    targets = ['Q13309']
#    for target in targets:
#        target_slice = df[df['Gene ID'].str.count(target)>0]
#        slice_vals = target_slice[names].values[0]
#        names_show = [ptna_underclip(name) for name in names]
#        error = [np.std(slice_vals[pd.Series(names_show)==name]) for name in names_show]
#        plt.ylabel('TMT Signal')
#        plt.xticks(rotation=52, ha="right", fontsize='small', rotation_mode="anchor")
#        plt.title(f"{target_slice['Gene Symbol'].values[0]}")
#        plt.bar(names_show, slice_vals, yerr=error, width=0.5, color='midnightblue')
#        plt.savefig(f'PDanalysis_results\\{folder}'+f'\\{target_slice["Gene Symbol"].values[0]}_sigsum.png',
#                    dpi=1200, bbox_inches='tight')
#        plt.clf()
#        plt.close()
#
#        plt.ylabel('TMT Signal')
#        plt.xticks(rotation=52, ha="right", fontsize='small', rotation_mode="anchor")
#        plt.title(f"{target_slice['Gene Symbol'].values[0]}")
#        plt.bar(names, slice_vals, width=0.5, color='midnightblue')
#        plt.savefig(f'PDanalysis_results\\{folder}'+f'\\{target_slice["Gene Symbol"].values[0]}_signal.png',
#                    dpi=1200, bbox_inches='tight')
#        plt.clf()
#        plt.close()


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

filepath = ROOT_DIR#[:-22] + '\\MS_Analyses_PD\\DSRD\\DSRD_4'
filename = '\\SKP2Global-(3).csv'

df = pd.read_csv(filepath+filename)

if os.path.exists('PDanalysis_results\\'+filename[:-4]):
    print(f'{filename[:-4]} exists')
else:
    os.mkdir('PDanalysis_results\\'+filename[:-4])

proteonalysis_main(filename, df)