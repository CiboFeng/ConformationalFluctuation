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
pywfile='pca.pyw'
npzfile='pca_distr.npz'
labels=['AA','HPS','Qch','Mid','Flx']
figname='pca'
nhi=50

"""Read Data and Calculate"""
perct=pyw(pywfile,'Percentage').reshape(-1)
pc=[]*len(labels)
for i in range(len(labels)):
    pc+=[pyw(pywfile,f'({labels[i]})')]

q=np.load(npzfile)
PC1=q['PC1']
P_PC1=q['P_PC1']
PC1_mean=q['PC1_mean']
PC1_std=q['PC1_std']
PC2=q['PC2']
P_PC2=q['P_PC2']
PC2_mean=q['PC2_mean']
PC2_std=q['PC2_std']
P_PC=q['P_PC']
# P_PC[np.isnan(P_PC)]=-np.inf

pc_lim=np.zeros((2,2))
pc_lim[0,0]=np.min(PC1)
pc_lim[0,1]=np.max(PC1)
pc_lim[1,0]=np.min(PC2)
pc_lim[1,1]=np.max(PC2)

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
    ax.plot(pc[i][0],pc[i][1],'.',linewidth=wth,alpha=0.2)
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
    ax.set_xlim(1.05*pc_lim[0,0]-0.05*pc_lim[0,1],1.05*pc_lim[0,1]-0.05*pc_lim[0,0])
    ax.set_ylim(1.05*pc_lim[1,0]-0.05*pc_lim[1,1],1.05*pc_lim[1,1]-0.05*pc_lim[1,0])
    if divmod(axord[i]-1,3)[0]==1:
        ax.set_xlabel(f'PC1 ({(perct[0]*100):.2f}%)',fontsize=size)
    if divmod(axord[i]-1,3)[1]==0:
        ax.set_ylabel(f'PC2 ({(perct[1]*100):.2f}%)',fontsize=size)
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
norm=colors.Normalize(vmin=0,vmax=np.max(P_PC))
ij=[(0,0),(0,1),(1,0),(1,1),(1,2)]

for i in range(len(labels)):
    ax=axs[ij[i]]
    levels=np.linspace(0.01*np.max(P_PC[i])+0.99*np.min(P_PC[i]),np.max(P_PC[i]),20)
    ax.contourf(PC1[i],PC2[i],P_PC[i].T,levels=levels,cmap=cmap)
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
    # ax.set_xlim(1.05*pc_lim[0,0]-0.05*pc_lim[0,1],1.05*pc_lim[0,1]-0.05*pc_lim[0,0])
    # ax.set_ylim(1.05*pc_lim[1,0]-0.05*pc_lim[1,1],1.05*pc_lim[1,1]-0.05*pc_lim[1,0])
    ax.set_xlim(-300,1300)
    ax.set_ylim(-500,500)
    if ij[i][0]==1:
        ax.set_xlabel(f'PC1 ({(perct[0]*100):.2f}%)',fontsize=size)
    if ij[i][1]==0:
        ax.set_ylabel(f'PC2 ({(perct[1]*100):.2f}%)',fontsize=size)
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

fig,axs=plt.subplots(2,2,figsize=(10,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
colors=['k',(0.5,0.5,0.5),'r','g','b']

ax=axs[0,0]
for i in range(len(labels)):
    ax.plot(PC1[i],P_PC1[i],color=colors[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('PC1',fontsize=size)
ax.set_ylabel('Probability Density',fontsize=size)
ax.set_xlim(-300,1300)
ax.legend(loc='upper right',fontsize=size)

ax=axs[0,1]
for i in range(len(labels)):
    ax.plot(PC2[i],P_PC2[i],color=colors[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('PC1',fontsize=size)
ax.set_ylabel('Probability Density',fontsize=size)
# ax.set_xlim(40,1500)

ax=axs[1,0]
ax.fill_between([0,len(labels)-1],[PC1_mean[0]-PC1_std[0],PC1_mean[0]-PC1_std[0]],[PC1_mean[0]+PC1_std[0],PC1_mean[0]+PC1_std[0]],color='k',linewidth=0,alpha=0.2)
for i in range(len(labels)):
    ax.errorbar(i,PC1_mean[i],yerr=PC1_std[i],
                    marker='.',markersize=0,markeredgewidth=wth,
                    ecolor=colors[i],elinewidth=wth,capsize=2*wth)
    ax.plot(i,PC1_mean[i],'.',color=colors[i],markersize=5*wth)
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
ax.set_xticks(np.arange(len(labels)),labels)
ax.set_ylabel('PC1',fontsize=size)
# ax.set_xlim(40,1500)

ax=axs[1,1]
ax.fill_between([0,len(labels)-1],[PC2_mean[0]-PC2_std[0],PC2_mean[0]-PC2_std[0]],[PC2_mean[0]+PC2_std[0],PC2_mean[0]+PC2_std[0]],color='k',linewidth=0,alpha=0.2)
for i in range(len(labels)):
    ax.errorbar(i,PC2_mean[i],yerr=PC2_std[i],
                    marker='.',markersize=0,markeredgewidth=wth,
                    ecolor=colors[i],elinewidth=wth,capsize=2*wth)
    ax.plot(i,PC2_mean[i],'.',color=colors[i],markersize=5*wth)
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
ax.set_xticks(np.arange(len(labels)),labels)
ax.set_ylabel('PC2',fontsize=size)
# ax.set_xlim(40,1500)

plt.tight_layout()
plt.savefig(f'{figname}_edge.png',format='png')
plt.savefig(f'{figname}_edge.pdf',format='pdf')
plt.show()
