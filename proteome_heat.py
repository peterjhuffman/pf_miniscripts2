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

savlink = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\37. INSM1 multidosing\exports\\'
x = ''
cmap_op = 'bwr'
save=True
spec = 2
vmn,vmx = -spec,spec

plt.rcParams.update({'font.size': 4})

filepath = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\results\\37. INSM1 multidosing\exports\\'
df = pd.read_excel(filepath+'PFEC_PH_20240624_ISDglobal_fx_prot_pure_analysisB_itx.xlsx')
#savlink = filepath +''

# df = df[df['type']!='crbn']
# df = df[df['type']!='AR targets']


print(df.shape)
df = df.dropna(subset=['type']).sort_values(by=['type'])
print(df.shape)
# df.to_csv(savlink+'heatmaps.csv',index=False)


heat_targets = ['d1 / ctrl Difference','d2 / ctrl Difference',
                'd3 / ctrl Difference','d4 / ctrl Difference',
                'd5 / ctrl Difference','d6 / ctrl Difference',
                'd7 / ctrl Difference','d8 / ctrl Difference']
#heat_targets = ['24hr_pf15 Difference','24hr_pf83 Difference']
#heat_targets = ['6hr_pf15 Difference','6hr_pf83 Difference']
heat_labels = [x[:-18] for x in heat_targets]



i=df[heat_targets].reset_index(drop=True).max().max()
j=df[heat_targets].reset_index(drop=True).min().min()
k = i/(i-j)


orig_cmap = cm.get_cmap(cmap_op)
# shifted_cmap = shiftedColorMap(orig_cmap, midpoint=1-k, name='shifted')

fig, ax = plt.subplots()
im = ax.imshow(df[heat_targets].T,cmap=cmap_op,interpolation='nearest', vmin=vmn,vmax=vmx)
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel('log2fc/DMSO', rotation=-90, va="bottom")


ax.set_xticks(np.arange(len(df['Gene Symbol'])))
ax.set_xticklabels(df['Gene Symbol'])
plt.xlabel('Gene Symbol')
ax.set_yticks(np.arange(len(heat_labels)))
ax.set_yticklabels(heat_labels)
# plt.ylabel('Treatment')
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

ax.set_title("Heatmap of target protein Log2FC/DMSO")
#ax.legend()
fig.tight_layout()

if save:
    plt.savefig(filepath+f'heat\\heatmap_full_{cmap_op}{x}{spec}.png', dpi=800)
plt.show()
plt.clf()


