"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir='jump'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
ntp=20
nat=163
nhi=200
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
lx=150.0
lxy=150*np.sqrt(2)
lz=500.0

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

d=np.zeros((nmd,nT,nat,nhi))
p_d=np.zeros((nmd,nT,nat,nhi))
dc=np.zeros((nmd,nT,nhi))
p_dc=np.zeros((nmd,nT,nhi))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Residue:')
        d[i,j]=q[::2]
        p_d[i,j]=q[1::2]
        q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Center')
        dc[i,j]=q[0]
        p_dc[i,j]=q[1]

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(20,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
seqidx=np.arange(nat)

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        for k in range(nat):
            rgb=cmap(norm(k))
            ax.plot(d[i,j,k],p_d[i,j,k],color=rgb,linewidth=wth)
        ax.plot(dc[i,j],p_dc[i,j],color='k',linewidth=wth)
        ax.plot([lx,lx],[0,max(np.max(p_d[i,j]),np.max(p_dc[i,j]))],color=(0.5,0.5,0.5),linestyle='--',linewidth=wth)
        ax.plot([lxy,lxy],[0,max(np.max(p_d[i,j]),np.max(p_dc[i,j]))],color=(0.5,0.5,0.5),linestyle='--',linewidth=wth)
        ax.plot([lz,lz],[0,max(np.max(p_d[i,j]),np.max(p_dc[i,j]))],color=(0.5,0.5,0.5),linestyle='--',linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        # ax.set_xlim(-7.5,157.5)
        if i==nmd-1:
            ax.set_xlabel('Displacement ($\\mathrm{\\AA}$)\nat %i K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Probability Density\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        ax.set_yscale('log')

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Residue Index',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()