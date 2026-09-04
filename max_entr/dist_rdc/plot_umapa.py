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
pywfile='umapa.pyw'
npzfile='umapa_distr.npz'
labels=['AA','HPS','Qch','Mid','Flx']
figname='umapa'
nhi=50

"""Read Data and Calculate"""
umap=[]*len(labels)
for i in range(len(labels)):
    umap+=[pyw(pywfile,f'({labels[i]})')]

q=np.load(npzfile)
UMAP1=q['UMAP1']
P_UMAP1=q['P_UMAP1']
UMAP2=q['UMAP2']
P_UMAP2=q['P_UMAP2']
P_UMAP=q['P_UMAP']
# P_UMAP[np.isnan(P_UMAP)]=-np.inf

umap_lim=np.zeros((2,2))
umap_lim[0,0]=np.min(UMAP1)
umap_lim[0,1]=np.max(UMAP1)
umap_lim[1,0]=np.min(UMAP2)
umap_lim[1,1]=np.max(UMAP2)

"""Plot"""
plt.figure(figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
axord=[1,2,4,5,6]

for i in range(len(labels)):
    ax=plt.subplot(2,3,axord[i])
    ax.plot(umap[i][0],umap[i][1],'.',linewidth=wth,alpha=0.2)
    ax.plot([],[],label=labels[i])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(1.05*umap_lim[0,0]-0.05*umap_lim[0,1],1.05*umap_lim[0,1]-0.05*umap_lim[0,0])
    ax.set_ylim(1.05*umap_lim[1,0]-0.05*umap_lim[1,1],1.05*umap_lim[1,1]-0.05*umap_lim[1,0])
    if divmod(axord[i]-1,3)[0]==1:
        ax.set_xlabel(f'UMAP1',fontsize=size)
    if divmod(axord[i]-1,3)[1]==0:
        ax.set_ylabel(f'UMAP2',fontsize=size)
    ax.legend(loc='best',fontsize=size,handlelength=0.0,handletextpad=0.0)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(2,3,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.85
cbarwth=0.03
xtick=0.1
ytick=20
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=np.max(P_UMAP))
ij=[(0,0),(0,1),(1,0),(1,1),(1,2)]

for i in range(len(labels)):
    ax=axs[ij[i]]
    levels=np.linspace(0.01*np.max(P_UMAP[i])+0.99*np.min(P_UMAP[i]),np.max(P_UMAP[i]),20)
    ax.contourf(UMAP1[i],UMAP2[i],P_UMAP[i].T,levels=levels,cmap=cmap)
    ax.plot([],[],label=labels[i])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(1.05*umap_lim[0,0]-0.05*umap_lim[0,1],1.05*umap_lim[0,1]-0.05*umap_lim[0,0])
    # ax.set_ylim(1.05*umap_lim[1,0]-0.05*umap_lim[1,1],1.05*umap_lim[1,1]-0.05*umap_lim[1,0])
    ax.set_xlim(-19,19)
    ax.set_ylim(-19,19)
    if ij[i][0]==1:
        ax.set_xlabel(f'UMAP1',fontsize=size)
    if ij[i][1]==0:
        ax.set_ylabel(f'UMAP2',fontsize=size)
    ax.legend(loc='best',fontsize=size,handlelength=0.0,handletextpad=0.0)
ax[0,2].set_axis_off()

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Density',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_distr.png',format='png')
plt.savefig(f'{figname}_distr.pdf',format='pdf')
plt.show()
