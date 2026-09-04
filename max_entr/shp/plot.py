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
# file='shp/1.37sgm_0.6sgm_itr20.npz'
# file='shp/1.5sgm_0.82sgm_itr20.npz'
file='shp/1.55sgm_1.0sgm_itr20.npz'
figname=file[:-4]
nhi=50

"""Read Data and Calculate"""
q=np.load(file)
Shp=q['Shp']
P1=q['P1']
P2=q['P2']

"""Plot"""
fig,axs=plt.subplots(6,6,figsize=(25,16))
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
norm=colors.Normalize(vmin=0,vmax=np.max(P2))
labels=['$R_1~(\\mathrm{\\AA})$','$R_2~(\\mathrm{\\AA})$','$R_3~(\\mathrm{\\AA})$','$R_\\mathrm{g}~(\\mathrm{\\AA})$','$\\mathit{\\Delta}$','$S$']

for i in range(6):
    for j in range(6):
        ax=axs[i,j]
        if i!=j:
            levels=np.linspace(0.01*np.max(P2[i,j])+0.99*np.min(P2[i,j]),np.max(P2[i,j]),20)
            ax.contourf(Shp[j],Shp[i],P2[j,i].T,levels=levels,cmap=cmap)
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
            ax.set_xlim(1.05*np.min(Shp[j])-0.05*np.max(Shp[j]),1.05*np.max(Shp[j])-0.05*np.min(Shp[j]))
            ax.set_ylim(1.05*np.min(Shp[i])-0.05*np.max(Shp[i]),1.05*np.max(Shp[i])-0.05*np.min(Shp[i]))
            if i==5:
                ax.set_xlabel(labels[j],fontsize=size)
            if j==0:
                ax.set_ylabel(labels[i],fontsize=size)
        else:
            ax.plot(Shp[i],P1[i],linewidth=wth)
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
            ax.set_xlim(1.05*np.min(Shp[i])-0.05*np.max(Shp[i]),1.05*np.max(Shp[i])-0.05*np.min(Shp[i]))
            # ax.set_ylim(1.05*pc_lim[1,0]-0.05*pc_lim[1,1],1.05*pc_lim[1,1]-0.05*pc_lim[1,0])
            if i==5:
                ax.set_xlabel(labels[j],fontsize=size)
            if j==0:
                ax.set_ylabel(labels[i],fontsize=size)

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
