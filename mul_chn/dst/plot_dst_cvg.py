"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
files=['dst_cvg/1.37sgm_0.6sgm.pyw',
       'dst_cvg/1.5sgm_0.82sgm.pyw',
       'dst_cvg/1.55sgm_1.0sgm.pyw']
figname='dst_cvg'
nhi=1000
nset=5
unit=0.0738e-3 ### 1g/(mol*AA)=0.0738mg/ml

"""Read Data and Calculate"""
z=np.zeros((len(files),nset,nhi))
dst=np.zeros((len(files),nset,nhi))
for i in range(len(files)):
    for j in range(nset):
        q=pyw(files[i],f'{j}-th')
        z[i,j]=q[0]
        dst[i,j]=q[1]

"""Plot"""
fig,ax=plt.subplots(3,1,figsize=(8,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nset-1)
labels=['Flex','Mid','Fix']

for i in range(len(files)):
    for j in range(nset):
        normc=norm(j)
        rgb=cmap(normc)
        ax[i].plot(z[i,j],dst[i,j]*unit,color=rgb,linewidth=wth)
    ax[i].autoscale()
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i].spines['bottom'].set_linewidth(wth)
    ax[i].spines['top'].set_linewidth(wth)
    ax[i].spines['left'].set_linewidth(wth)
    ax[i].spines['right'].set_linewidth(wth)
    # ax[i].set_xlim(-4,84)
    if i==len(files)-1:
        ax[i].set_xlabel('z ($\\mathrm{\\AA}$)',fontsize=size)
    ax[i].set_ylabel('Density (g/ml)\nof %s Model'%labels[i],fontsize=size)

fig.subplots_adjust(right=0.85)
plt.tight_layout(rect=[0,0,0.85,1])
top=ax[0].get_position().y1
bot=ax[2].get_position().y0
cbar_ax=fig.add_axes([0.88,bot,0.03,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Set',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.set_ticks([0,1,2,3,4])
cbar.set_ticklabels([1,2,3,4,5])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()