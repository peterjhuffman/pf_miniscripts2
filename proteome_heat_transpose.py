# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:11:44 2023

@author: HUFFMP
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

savlink = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\AHI insm1 rtsMS3 global\\exports\\'
x = '_STAT_'
cmap_op = 'bwr'
save=True
spec = 3
vmn,vmx = -spec,spec

plt.rcParams.update({'font.size': 4})

filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\AHI insm1 rtsMS3 global\\exports\\'
df = pd.read_excel(filepath+'PFAS_PH_20240405_AHIglobal_rtsMS3_prot_analysis_itxome.xlsx')
#savlink = filepath +''

# df = df[df['type']!='crbn']
# df = df[df['type']!='AR targets']


print(df.shape)
df = df.dropna(subset=['type2']).sort_values(by=['summative'], ascending=True)
print(df.shape)
# df.to_csv(savlink+'heatmaps.csv',index=False)


# heat_targets = ['4h Low / 4h DMSO Difference','4h High / 4h DMSO Difference',
#                 '18h Low / 18h DMSO Difference','18h High / 18h DMSO Difference']
heat_targets = ['4 lo / 4 dmso Difference','4 hi / 4 dmso Difference',
                '18 lo / 18 dmso Difference','18 hi / 18 dmso Difference']

heat_labels = [x.split(' ')[0]+x.split(' ')[1] for x in heat_targets]



i=df[heat_targets].reset_index(drop=True).max().max()
j=df[heat_targets].reset_index(drop=True).min().min()
k = i/(i-j)


orig_cmap = cm.get_cmap(cmap_op)
# shifted_cmap = shiftedColorMap(orig_cmap, midpoint=1-k, name='shifted')

fig, ax = plt.subplots()
im = ax.imshow(df[heat_targets].T,cmap=cmap_op,interpolation='nearest', vmin=vmn,vmax=vmx)
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_xlabel('log2fc/DMSO',  va="top")


ax.set_xticks(np.arange(len(df['Gene Symbol'])))
ax.set_xticklabels(df['Gene Symbol'], rotation=90)
ax.set_yticks(np.arange(len(heat_labels)))
ax.set_yticklabels(heat_labels)
plt.setp(ax.get_yticklabels(), ha="right",
         rotation_mode="anchor")

ax.set_title("Heatmap of target protein Log2FC/DMSO")
#ax.legend()
fig.tight_layout()

if save:
    plt.savefig(filepath+f'heat\\heatmap_full_{cmap_op}{x}{spec}.png', dpi=800)
plt.show()
plt.clf()


