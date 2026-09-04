"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
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
dir='msd'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nfr=10000
nat=163
cut=20

"""Read Data"""
t=np.zeros((nmd,nfr,nat))
d=np.zeros((nmd,nfr,nat))
for i in range(nmd):
    q=pyw(f'{dir}/{mdls[i]}_dyn.pyw','Mean-squared')
    for j in range(nfr):
        for k in range(nat):
            t[i,j,k]=q[0,j*nat+k]
            d[i,j,k]=q[2,j*nat+k]
dm=np.mean(d,axis=-1)
dmax=np.mean(d[:,200:5000],axis=1)
dmid=np.sum(d[:,:10]*np.arange(nfr-1,nfr-11,-1)[np.newaxis,:,np.newaxis],axis=1)/np.sum(np.arange(nfr-1,nfr-11,-1))

"""Fit"""
def func(x,a,b):
    return a*x**b

Par=np.ones((nmd,2))
for i in range(nmd):
    par,cov=curve_fit(func,t[i,:cut,0],dm[i,:cut])
    Par[i]=par

"""Plot"""
fig,axs=plt.subplots(2,3,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.82
cbarwth=0.03
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
clrs=['r','g','b']

for i in range(nmd):
    ax=axs[0,i]
    for j in range(nat):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(t[i,:,j],d[i,:,j],color=rgb,linewidth=wth,alpha=0.5)
    ax.plot(t[i,:,0],dm[i,:],color='k',linewidth=wth)
    if i==0:
        ax.plot([],[],color='k',linewidth=wth,label='Mean')
    ax.plot(t[i,1:cut,0],3*Par[i,0]*t[i,1:cut,0]**Par[i,1],color=(0.5,0.5,0.5),linestyle=':',linewidth=wth,label='$\\langle{}\\mathrm{\\Delta}r^2\\rangle{}\\propto{}\\mathrm{\\Delta}t^{%.2f}$'%Par[i,1])
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
    ax.set_xlabel('Lag Time (ns)',fontsize=size)
    ax.set_ylabel('MSD of %s Model ($\\mathrm{\\AA}^2$)'%labels[i],fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(40,1500)
    ax.legend(loc='lower right',fontsize=size)

ax=axs[1,0]
for i in range(nmd):
    ax.plot(t[i,:,0],dm[i,:],color=clrs[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('Lag Time (ns)',fontsize=size)
ax.set_ylabel('MSD of %s Model ($\\mathrm{\\AA}^2$)'%labels[i],fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(40,1500)
ax.legend(loc='lower right',fontsize=size)

ax=axs[1,1]
for i in range(nmd):
    ax.plot(np.arange(nat),dmid[i],color=clrs[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('Residue Index',fontsize=size)
ax.set_ylabel('Middle MSD of %s Model ($\\mathrm{\\AA}^2$)',fontsize=size)
ax.legend(loc='upper center',fontsize=size)

ax=axs[1,2]
for i in range(nmd):
    ax.plot(np.arange(nat),dmax[i],color=clrs[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('Residue Index',fontsize=size)
ax.set_ylabel('Maximum MSD of %s Model ($\\mathrm{\\AA}^2$)',fontsize=size)
ax.legend(loc='upper center',fontsize=size)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Residue Index',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()