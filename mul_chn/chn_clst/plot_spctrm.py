"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
from scipy.optimize import curve_fit
from scipy import optimize
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
file='proc.npz'
figname='chn_clst_spctrm'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nrep=20
nfr=10000
nch=50

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

q=np.load(file)
t=q['t']
nclst=q['nclst']
nclst_mean=q['nclst_mean']
nclst_std=q['nclst_std']
treq=q['treq']
treq_mean=q['treq_mean']
treq_std=q['treq_std']
nisol=q['nisol']
nisol_mean=q['nisol_mean']
nisol_std=q['nisol_std']
nsq=q['nsq']
dnsq=q['dnsq']
dnsq_mean=q['dnsq_mean']
dnsq_std=q['dnsq_std']
dnmin=q['dnmin']
dnmin_mean=q['dnmin_mean']
dnmin_std=q['dnmin_std']
ndnmin=np.sum(dnmin!=0,axis=-1)

nclst_means=nclst_mean[:,0:1]
nclst_meane=nclst_mean[:,-1:]
nclst_mean=(nclst_mean-nclst_means)/(nclst_meane-nclst_means)

"""Fit nclst"""
# def exp_expa(x,f,k,alf=1e-3,no_neg=False):
#     """Small alf: good fit, oscillative solvation, sensitive to noise. Large alf: smooth solvation"""
#     from scipy.optimize import lsq_linear
#     nk=len(k)

#     w=np.zeros(nk)
#     w[1:-1]=(k[2:]-k[:-2])/2
#     w[0]=(k[1]-k[0])/2
#     w[-1]=(k[-1]-k[-2])/2
#     E=np.exp(-x[:,None]*k[None,:])
#     if alf>0:
#         reg=np.eye(nk)*np.sqrt(alf)
#         A=np.vstack([E,reg])
#         b=np.concatenate([f,np.zeros(nk)])
#     else:
#         A=E
#         b=f
    
#     if no_neg:
#         res=lsq_linear(A,b,bounds=(0,np.inf),method='trf',tol=1e-10)
#     else:
#         res=lsq_linear(A,b,bounds=(-np.inf,np.inf),method='trf',tol=1e-10)
#     m_sol=res.x
#     coeff=m_sol/w
#     f_rec=np.dot(E,m_sol)
    
#     return coeff,f_rec

def exp_expa(x,f,k,alf=1e-3,no_neg=False):
    """Small alf: good fit, oscillative solvation, sensitive to noise. Large alf: smooth solvation"""
    from scipy.optimize import lsq_linear
    nk=len(k)

    E=np.exp(-x[:,None]*k[None,:])
    if alf>0:
        reg=np.eye(nk)*np.sqrt(alf)
        A=np.vstack([E,reg])
        b=np.concatenate([f,np.zeros(nk)])
    else:
        A=E
        b=f
    
    if no_neg:
        res=lsq_linear(A,b,bounds=(0,np.inf),method='trf',tol=1e-10)
    else:
        res=lsq_linear(A,b,bounds=(-np.inf,np.inf),method='trf',tol=1e-10)
    coeff=res.x
    f_rec=np.dot(E,coeff)
    
    return coeff,f_rec

nk=100
kintra=np.logspace(-2,1,nk)
c=np.zeros((nmd,nk))
f=np.zeros((nmd,nfr))
for i in range(nmd):
    c[i],f[i]=exp_expa(t,1-nclst_mean[i],kintra,1e-3,no_neg=True)

"""Plot"""
fig,axs=plt.subplots(1,2,figsize=(12,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.85
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=axs[0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t,nclst_mean[i],color=rgb,linewidth=wth)
    ax.plot(t,1-f[i],color='k',linestyle='--',linewidth=wth)
ax.plot([],[],'k--',linewidth=wth,label='Reproduced')
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
# ax.set_xlim(-4,84)
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Intra-chain\nContact Formation',fontsize=size)
ax.legend(loc='best',fontsize=size)

ax=axs[1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(kintra[1:],c[i,1:],color=rgb,linewidth=wth,)
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
# ax.set_ylim(-10,100)
ax.set_xlabel('Rate (ns$^{-1}$)',fontsize=size)
ax.set_ylabel('Amplitude (ns)',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
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
cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()