"""Import Modules"""
from scipy.optimize import curve_fit
from scipy import optimize
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
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Fix','Mid','Flex']
nmd=len(mdls)
file='cont_frac_kin_annl.pyw'
figname=file[:-4]
nfr=10000

"""Read Data and Calculate"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

t=np.zeros((nmd,nfr))
Qintra=np.zeros((nmd,nfr))
Qinter=np.zeros((nmd,nfr))
for i in range(nmd):
    q=pyw(file,f'({Mdls[i]})')
    t[i]=q[0]
    Qintra[i]=q[1]
    Qinter[i]=q[2]

# t=t[:,1:]
# Qintra=Qintra[:,1:]
# Qinter=Qinter[:,1:]
Qintras=Qintra[:,0:1]
Qintrae=np.mean(Qintra[:,-1000:],axis=-1,keepdims=True)
Qintra=(Qintra-Qintras)/(Qintrae-Qintras)
Qinters=Qinter[:,0:1]
Qintere=np.mean(Qinter[:,-1000:],axis=-1,keepdims=True)
Qinter=(Qinter-Qinters)/(Qintere-Qinters)

""""""
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
cintra=np.zeros((nmd,nk))
fintra=np.zeros((nmd,nfr))
kinter=np.logspace(-2,1,nk)
cinter=np.zeros((nmd,nk))
finter=np.zeros((nmd,nfr))
for i in range(nmd):
    cintra[i],fintra[i]=exp_expa(t[i],1-Qintra[i],kintra,1e-3,no_neg=True)
    cinter[i],finter[i]=exp_expa(t[i],1-Qinter[i],kinter,1e-3,no_neg=True)

"""Plot"""
fig,axs=plt.subplots(2,2,figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.83
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=axs[0,0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i],Qintra[i],color=rgb,linewidth=wth)
    ax.plot(t[i],1-fintra[i],color='k',linestyle='--',linewidth=wth)
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

ax=axs[0,1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i],Qinter[i],color=rgb,linewidth=wth)
    ax.plot(t[i],1-finter[i],color='k',linestyle='--',linewidth=wth)
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
# ax.set_yscale('symlog')
# ax.set_ylim(-10,0)
# ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Inter-chain\nContact Formation',fontsize=size)

ax=axs[1,0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(kintra[1:],cintra[i,1:],color=rgb,linewidth=wth)
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

ax=axs[1,1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(kinter[1:],cinter[i,1:],color=rgb,linewidth=wth)
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
# ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
ax.set_xlabel('Rate (ns$^{-1}$)',fontsize=size)
ax.set_ylabel('Amplitude (ns)',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

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
cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_spctrm.png',format='png')
plt.savefig(f'{figname}_spctrm.pdf',format='pdf')
plt.show()

# fig=plt.figure(figsize=(10,8))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# clrs=['r','g','b']

# ax=plt.subplot(2,2,1)
# for i in range(nmd):
#     ax.plot(t[i],Qintra[i],color=clrs[i],linewidth=wth,label=mdls[i])
#     ax.plot(t[i],1-fintra[i],color='k',linestyle='--',linewidth=wth)
# ax.plot([],[],'k--',linewidth=wth,label='Reproduced')
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_xlim(-4,84)
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Intra-chain\nContact Formation',fontsize=size)
# ax.legend(loc='best',fontsize=size)

# ax=plt.subplot(2,2,2)
# for i in range(nmd):
#     ax.plot(t[i],Qinter[i],color=clrs[i],linewidth=wth,label=mdls[i])
#     ax.plot(t[i],1-finter[i],color='k',linestyle='--',linewidth=wth)
# ax.plot([],[],'k--',linewidth=wth,label='Reproduced')
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_yscale('symlog')
# # ax.set_ylim(-10,0)
# # ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Inter-chain\nContact Formation',fontsize=size)

# ax=plt.subplot(2,2,3)
# for i in range(nmd):
#     ax.plot(kintra[1:],cintra[i,1:],color=clrs[i],linewidth=wth,)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_ylim(-10,100)
# ax.set_xlabel('Rate (ns$^{-1}$)',fontsize=size)
# ax.set_ylabel('Amplitude (ns)',fontsize=size)
# # ax.legend(loc='best',fontsize=0.5*size)

# ax=plt.subplot(2,2,4)
# for i in range(nmd):
#     ax.plot(kinter[1:],cinter[i,1:],color=clrs[i],linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_ylim(-10,100)
# # ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
# ax.set_xlabel('Rate (ns$^{-1}$)',fontsize=size)
# ax.set_ylabel('Amplitude (ns)',fontsize=size)
# # ax.legend(loc='best',fontsize=0.5*size)

# plt.tight_layout()
# plt.savefig(f'{figname}_spctrm.png',format='png')
# plt.savefig(f'{figname}_spctrm.pdf',format='pdf')
# plt.show()