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
mdls=['rg_kin/1.37sgm_0.6sgm',
      'rg_kin/1.5sgm_0.82sgm',
      'rg_kin/1.55sgm_1.0sgm']
figname='rg_kin'
nrep=20
nch=50
dt=10e-6*10000 ### ns

"""Read Data and Calculate"""
rg=[]
for i in range(len(mdls)):
    rg.append([])
    for j in range(nrep):
        rg[-1].append(pyw(f'{mdls[i]}_{j}.pyw','Radius')[:,:-1])

nfr=np.zeros((len(mdls),nrep),dtype=int)
for i in range(len(mdls)):
     for j in range(nrep):
          nfr[i,j]=np.shape(rg[i][j])[1]
nfrmax=np.max(nfr)
t=np.arange(nfrmax)*dt

rg_mean=np.zeros((len(mdls),nfrmax))
rg_std=np.zeros((len(mdls),nfrmax))
for i in range(len(mdls)):
    for j in range(nfrmax):
        rg_reshp=[]
        for k in range(nrep):
            try:
                rg_reshp+=rg[i][k][:,j].tolist()
            except Exception as err:
                pass
        rg_mean[i,j]=np.mean(rg_reshp)
        rg_std[i,j]=np.std(rg_reshp)

"""Plot"""
fig=plt.figure(figsize=(12,8))
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
norm=colors.Normalize(vmin=0,vmax=nrep*nch)
clrs=[(0.82,0.66,0.82),(0.56,0.90,0.94),(0.63,0.76,0.58)]
lbls=['Flex','Mid','Fix']
clrs_=['r','b']

for i in range(len(mdls)):
    ax=plt.subplot(2,2,i+1)
    for j in range(nrep):
        for k in range(nch):
            normc=norm(j*nch+k)
            rgb=cmap(normc)
            ax.plot(np.arange(len(rg[i][j][k]))*dt,rg[i][j][k],color=rgb,linewidth=wth,alpha=0.5)
    ax.plot(t,rg_mean[i],color='k',linewidth=wth,label='Mean')
    # ax.plot(t,nch*np.exp(-Par[i]*t),color='k',linestyle=':',linewidth=wth,label='$y=50\\mathrm{e}^{-%.2fx}$'%Par[i])
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
    ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_ylim(-2.5,52.5)
    if divmod(i,2)[0]==1:
        ax.set_xlabel('Time (ns)',fontsize=size)
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)\n of %s Model'%lbls[i],fontsize=size)
    ax.legend(loc='upper left',fontsize=size)

ax=plt.subplot(2,2,4)
ax2=ax.twinx()
for i in range(len(mdls)):
    ax.plot(t,rg_mean[i],color=clrs_[0],linewidth=3*wth)
    ax.plot(t,rg_mean[i],color=clrs[i],linewidth=2*wth,label=lbls[i])
    ax2.plot(t,rg_std[i],color=clrs_[1],linewidth=3*wth)
    ax2.plot(t,rg_std[i],color=clrs[i],linewidth=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=clrs_[0],labelcolor=clrs_[0])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=clrs_[0],labelcolor=clrs_[0])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(0)
ax.spines['left'].set_color(clrs_[0])
ax.yaxis.label.set_color(clrs_[0]) 
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim(-0.72,28.32)
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Radius of Gyration\nMean ($\\mathrm{\\AA}$)',fontsize=size)
ax.legend(loc='upper left',fontsize=size)
ax2.minorticks_on()
ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=clrs_[1],labelcolor=clrs_[1])
ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=clrs_[1],labelcolor=clrs_[1])
# ax2.xaxis.set_major_locator(MultipleLocator(xtick))
# ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.yaxis.set_major_locator(MultipleLocator(ytick))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(wth)
ax2.spines['right'].set_color(clrs_[1])
ax2.yaxis.label.set_color(clrs_[1])
ax2.set_ylim(-0.72,28.32)
ax2.set_ylabel('Radius of Gyration\nDeviation ($\\mathrm{\\AA}$)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()