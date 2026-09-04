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
npzfile='pca_shp.npz'
figname='pca_shp'
labels=['AA','HPS','Qch','Mid','Flx']
nmd=len(labels)-2
show='Mid'

"""Read Data and Calculate"""
data=np.load(npzfile)
cor=data['cor']
q_mean=data['q_mean']
q_std=data['q_std']
Q=data['Q']
P1=data['P1']
P2=data['P2']
nq=np.shape(cor)[1]

"""Plot"""
plt.figure(figsize=(30,25))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=np.max(P2))
axlabels=['PC1','PC2','$R_1~(\\mathrm{\\AA})$','$R_2~(\\mathrm{\\AA})$','$R_3~(\\mathrm{\\AA})$','$R_\\mathrm{g}~(\\mathrm{\\AA})$','$\\mathit{\\Delta}$','$S$']

i=labels.index(show)
for j in range(nq):
    for k in range(nq):
        ax=plt.subplot(nq,nq,j*nq+k+1)
        if j!=k:
            levels=np.linspace(0.01*np.max(P2[i,j,k])+0.99*np.min(P2[i,j,k]),np.max(P2[i,j,k]),20)
            ax.contourf(Q[i,k],Q[i,j],P2[i,k,j].T,levels=levels,cmap=cmap)
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
            if j==nq-1:
                ax.set_xlabel(axlabels[k],fontsize=size)
            if k==0:
                ax.set_ylabel(axlabels[j],fontsize=size)
        else:
            ax.plot(Q[i,j],P1[i,j],linewidth=wth)
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
            ax.set_xlim(1.05*np.min(Q[i,j])-0.05*np.max(Q[i,j]),1.05*np.max(Q[i,j])-0.05*np.min(Q[i,j]))
            # ax.set_ylim(1.05*pc_lim[1,0]-0.05*pc_lim[1,1],1.05*pc_lim[1,1]-0.05*pc_lim[1,0])
            if j==nq-1:
                ax.set_xlabel(labels[k],fontsize=size)
            if k==0:
                ax.set_ylabel(labels[j],fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_{show.lower()}.png',format='png')
plt.savefig(f'{figname}_{show.lower()}.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(2,3,figsize=(15,9))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.85
cbarwth=0.03
xtick=0.1
ytick=20
axord=[1,2,4,5,6]
color=[(0,0,1),(1,1,1),(1,0,0)]
nodes=[0/2,1/2,2/2]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=-1,vmax=1)
axlabels=['PC1','PC2','$R_1$','$R_2$','$R_3$','$R_\\mathrm{g}$','$\\mathit{\\Delta}$','$S$']
ij=[(0,0),(0,1),(1,0),(1,1),(1,2)]

for i in range(len(labels)):
    ax=axs[ij[i]]
    ax.imshow(cor[i],cmap=cmap,norm=norm)
    ax.invert_yaxis()
    ax.plot([],[],label=labels[i])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    # ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    if ij[i][0]==1:
        ax.set_xticks(np.arange(nq),axlabels,rotation=45)
    else:
        ax.set_xticks(np.arange(nq),[])
    if ij[i][1]==0:
        ax.set_yticks(np.arange(nq),axlabels)
    else:
        ax.set_yticks(np.arange(nq),[])
    ax.legend(loc='upper left',fontsize=size,handlelength=0.0,handletextpad=0)
ax[0,2].set_axis_off()

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Correlation Coefficient',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()