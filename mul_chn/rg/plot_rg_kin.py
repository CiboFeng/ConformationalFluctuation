"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir='rg_kin_annl'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
figname=dir
nfr=10000
nrep=20
nch=50
dt=10e-6*10000 ### ns

"""Read Data and Calculate"""
rg0_mean=np.zeros(nmd)
rg0_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg0_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg0_std[i]=q[0,0]
rg0_flu=np.round(np.log10(rg0_std/rg0_mean),1)

t=np.arange(nfr)*dt
rg=np.zeros((nmd,nrep,nfr,nch))
for i in range(nmd):
    for j in range(nrep):
        rg[i,j]=pyw(f'{dir}/{mdls[i]}_{j}.pyw','Radius').T
rg_mean=np.mean(rg,axis=(1,3))
rg_std=np.std(rg,axis=(1,3))
print(rg_mean[:,:2])
print(rg_std[:,:2])

"""Plot"""
fig,axs=plt.subplots(2,2,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
cmap_=plt.cm.rainbow
norm_=colors.Normalize(vmin=0,vmax=nrep*nch)
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
clrs_=['r','g']

for i in range(nmd):
    ax=axs[divmod(i,2)]
    for j in range(nrep):
        for k in range(1):
            rgb=cmap_(norm_(j*nch+k))
            ax.plot(t,rg[i,j,:,k],color=rgb,linewidth=wth,alpha=0.5)
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
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)\n of Flu.=%.2f'%10**rg0_flu[i],fontsize=size)
    ax.legend(loc='upper left',fontsize=size)

ax=axs[1,1]
ax2=ax.twinx()
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t,rg_mean[i],linestyle='-',color=clrs_[0],linewidth=3*wth)
    ax.plot(t,rg_mean[i],linestyle='-',color=rgb,linewidth=2*wth)
    ax.plot(t,rg_std[i],linestyle=':',color=clrs_[0],linewidth=3*wth)
    ax.plot(t,rg_std[i],linestyle='-',color=rgb,linewidth=2*wth)
    ax2.plot(t,rg_std[i]/rg_mean[i],color=clrs_[1],linewidth=3*wth)
    ax2.plot(t,rg_std[i]/rg_mean[i],color=rgb,linewidth=2*wth)
    # ax.plot(t[:100],rg_mean[i,:100],linestyle='-',color=clrs_[0],linewidth=3*wth)
    # ax.plot(t[:100],rg_mean[i,:100],linestyle='-',color=rgb,linewidth=2*wth)
    # ax.plot(t[:100],rg_std[i,:100],linestyle=':',color=clrs_[0],linewidth=3*wth)
    # ax.plot(t[:100],rg_std[i,:100],linestyle='-',color=rgb,linewidth=2*wth)
    # ax2.plot(t[:100],rg_std[i,:100]/rg_mean[i,:100],color=clrs_[1],linewidth=3*wth)
    # ax2.plot(t[:100],rg_std[i,:100]/rg_mean[i,:100],color=rgb,linewidth=2*wth)
ax2.plot([],[],linestyle='-',color=clrs_[0],linewidth=wth,label='Mean')
ax2.plot([],[],linestyle=':',color=clrs_[0],linewidth=wth,label='Deviation')
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
# ax.set_ylim(-0.72,28.32)
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Radius of Gyration\nMean & Deviation ($\\mathrm{\\AA}$)',fontsize=size)
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
ax2.set_ylim(-0.02,0.36)
ax2.set_ylabel('Radius of Gyration\nFluctuation',fontsize=size)
ax2.legend(loc='upper left',fontsize=0.8*size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Structure Fluctuation',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
cbar.set_ticklabels([f'{10**rg0_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()