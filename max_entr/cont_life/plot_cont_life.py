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
dir='cont_life'
# dir='cont_life_lin'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nat=163
nhi=20

"""Read Data and Calculate"""
T0=np.zeros((nmd,nat,nat,nhi))
p_T0=np.zeros((nmd,nat,nat,nhi))
T0_mean=np.zeros((nmd,nat,nat))
T0_std=np.zeros((nmd,nat,nat))
T1=np.zeros((nmd,nat,nat,nhi))
p_T1=np.zeros((nmd,nat,nat,nhi))
T1_mean=np.zeros((nmd,nat,nat))
T1_std=np.zeros((nmd,nat,nat))
for i in range(nmd):
    q=np.load(f'{dir}/{mdls[i]}_dyn.npz')
    T0[i]=q['T0']
    p_T0[i]=q['p_T0']
    T0_mean[i]=q['T0_mean']
    T0_std[i]=q['T0_std']
    T1[i]=q['T1']
    p_T1[i]=q['p_T1']
    T1_mean[i]=q['T1_mean']
    T1_std[i]=q['T1_std']

"""Plot"""
fig,axs=plt.subplots(3,2,figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
label0=['Qch','Mid','Flx']
npair=10
Dseq=np.arange(nat-1,0,-round(nat/npair))
IJ=[]
for i in range(npair):
    I=random.randint(0,nat-1-Dseq[i])
    J=I+Dseq[i]
    IJ+=[(I,J)]

for i in range(nmd):
    ax=axs[i,0]
    for l,m in IJ:
        normc=norm(abs(l-m))
        rgb=cmap(normc)
        ax.plot(T0[i,l,m],p_T0[i,l,m],color=rgb,linewidth=wth)
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
    # ax.set_xscale('log')
    ax.set_xticks([-1,0,1,2,3],['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
    if i==nmd-1:
        ax.set_xlabel('Life Time of 0 (ns)',fontsize=size)
    ax.set_ylabel('Probability Density\nof %s Model ($\\mathrm{\\AA}^{-1}$)'%label0[i],fontsize=size)

    ax=axs[i,1]
    for l,m in IJ:
        normc=norm(abs(l-m))
        rgb=cmap(normc)
        ax.plot(T1[i,l,m],p_T1[i,l,m],color=rgb,linewidth=wth)
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
    # ax.set_xscale('log')
    if i==nmd-1:
        ax.set_xlabel('Life Time of 1 (ns)',fontsize=size)
    # ax.set_ylabel('Probability Density\nof %s Model ($\\mathrm{\\AA}^{-1}$)'%label0[i],fontsize=size)
    ax.set_xticks([-1,0,1,2,3],['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()