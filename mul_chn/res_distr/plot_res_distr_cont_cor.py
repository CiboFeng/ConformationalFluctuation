"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir1='../dst/dst'
dir='res_distr'
dir2='../cont/cont_slab'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
ntp=20
nat=163
nhi=250
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

# """Read Data"""
# rg_mean=np.zeros(nmd)
# rg_std=np.zeros(nmd)
# for i in range(nmd):
#     q=pyw(f'{dir0}/{mdls[i]}_{Ts[0]}.pyw','Mean')
#     rg_mean[i]=q[0,0]
#     q=pyw(f'{dir0}/{mdls[i]}_{Ts[0]}.pyw','Deviation')
#     rg_std[i]=q[0,0]
# rg_flu=np.round(np.log10(rg_std/rg_mean),1)

# z_seq=np.zeros((nmd,nT,nat,nhi))
# pz_seq=np.zeros((nmd,nT,nat,nhi))
# zabs_seq_mean=np.zeros((nmd,nT,nat))
# zabs_seq_std=np.zeros((nmd,nT,nat))
# z_tp=np.zeros((nmd,nT,ntp,nhi))
# pz_tp=np.zeros((nmd,nT,ntp,nhi))
# zabs_tp_mean=np.zeros((nmd,nT,ntp))
# zabs_tp_std=np.zeros((nmd,nT,ntp))
# for i in range(nmd):
#     for j in range(nT):
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Residue Index,')
#         z_seq[i,j]=q[1].reshape(nat,nhi)
#         pz_seq[i,j]=q[2].reshape(nat,nhi)
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Mean of Radial Position of Each Residue:')
#         zabs_seq_mean[i,j]=q[0]
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Standard Deviation of Radial Position of Each Residue:')
#         zabs_seq_std[i,j]=q[0]
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Residue Type,')
#         z_tp[i,j]=q[1].reshape(ntp,nhi)
#         pz_tp[i,j]=q[2].reshape(ntp,nhi)
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Mean of Radial Position of Each Residue Type:')
#         zabs_tp_mean[i,j]=q[0]
#         q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Standard Deviation of Radial Position of Each Residue Type:')
#         zabs_tp_std[i,j]=q[0]

# """Renormalize"""
# p_seq=pz_seq/np.sum(pz_seq,axis=2,keepdims=True)
# p_seq/=np.max(p_seq,axis=2,keepdims=True)
# popu=np.zeros(ntp,dtype=int)
# for i in range(nat):
#     popu[seq[i]]+=1
# # p_tp=p_r_tp*popu[np.newaxis,np.newaxis,:,np.newaxis]
# p_tp=pz_tp.copy()
# p_tp/=np.sum(p_tp,axis=2,keepdims=True)
# p_tp/=np.max(p_tp,axis=2,keepdims=True)

# """Sort"""
# idx_sort=[res.index(i) for i in res_sort]
# z_tp=z_tp[:,:,idx_sort]
# pz_tp=pz_tp[:,:,idx_sort]
# p_tp=p_tp[:,:,idx_sort]
# zabs_tp_mean=zabs_tp_mean[:,:,idx_sort]
# zabs_tp_std=zabs_tp_std[:,:,idx_sort]

# zabs_tp_mean[zabs_tp_mean==0.0]=np.nan
# zabs_tp_std[zabs_tp_std==0.0]=np.nan

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]

z_seq=np.zeros((nmd,nT,nhi))
pz_seq=np.zeros((nmd,nT,nat,nhi))
zabs_seq_mean=np.zeros((nmd,nT,nat))
zabs_seq_std=np.zeros((nmd,nT,nat))
z_tp=np.zeros((nmd,nT,nhi))
pz_tp=np.zeros((nmd,nT,ntp,nhi))
zabs_tp_mean=np.zeros((nmd,nT,ntp))
zabs_tp_std=np.zeros((nmd,nT,ntp))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        z_seq[i,j]=q['Z']
        pz_seq[i,j]=q['pZ'].T
        zabs_seq_mean[i,j]=q['zabs_mean']
        zabs_seq_std[i,j]=q['zabs_std']
        z_tp[i,j]=q['Z_tp']
        pz_tp[i,j]=q['pZ_tp'].T
        zabs_tp_mean[i,j]=q['zabs_tp_mean']
        zabs_tp_std[i,j]=q['zabs_tp_std']

"""Renormalize"""
pz_seq=pz_seq/np.sum(pz_seq,axis=2,keepdims=True)
pz_seq/=np.max(pz_seq,axis=2,keepdims=True)
pz_tp=pz_tp/np.sum(pz_tp,axis=2,keepdims=True)
popu=np.zeros(ntp,dtype=int)
for i in range(nat):
    popu[seq[i]]+=1
for i in range(ntp):
    if popu[i]>0:
        pz_tp[:,:,i]/=popu[np.newaxis,np.newaxis,i,np.newaxis]
pz_tp/=np.max(pz_tp,axis=2,keepdims=True)

"""Sort"""
idx_sort=[res.index(i) for i in res_sort]
pz_tp=pz_tp[:,:,idx_sort]
zabs_tp_mean=zabs_tp_mean[:,:,idx_sort]
zabs_tp_std=zabs_tp_std[:,:,idx_sort]

zabs_tp_mean[zabs_tp_mean==0.0]=np.nan
zabs_tp_std[zabs_tp_std==0.0]=np.nan

"""Read Contact Data"""
P_intra=np.zeros((nmd,nT,nat,nat))
P_inter=np.zeros((nmd,nT,nat,nat))
p_intra=np.zeros((nmd,nT,ntp,ntp))
p_inter=np.zeros((nmd,nT,ntp,ntp))
num=np.zeros((ntp,ntp),dtype=int)
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir2}/{mdls[i]}_{Ts[j]}.npz')
        P_intra[i,j]=q['Pintra']
        P_inter[i,j]=q['Pinter']
for i in range(4):
    for j in range(nat-i):
        P_intra[:,:,j,i+j]=0.0
for i in range(nat):
    for j in range(i+1):
        p_intra[:,:,seq[i],seq[j]]+=P_intra[:,:,i,j]
        p_inter[:,:,seq[i],seq[j]]+=P_inter[:,:,i,j]
        num[seq[i],seq[j]]+=1
        if seq[i]!=seq[j]:
            p_intra[:,:,seq[j],seq[i]]+=P_intra[:,:,i,j]
            p_inter[:,:,seq[j],seq[i]]+=P_inter[:,:,i,j]
            num[seq[j],seq[i]]+=1
num[num==0]=1
# num=np.ones((nres,nres),dtype=int)
p_intra=np.sum(p_intra/num[np.newaxis,np.newaxis],axis=2)[:,:,idx_sort]
p_inter=np.sum(p_inter/num[np.newaxis,np.newaxis],axis=2)[:,:,idx_sort]

"""Plot"""
fig,axs=plt.subplots(2,nT,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
typeidx=np.arange(ntp)

for i in range(nT):
    ax=axs[0,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.scatter(zabs_tp_mean[j,i],p_intra[j,i],color=rgb,s=20*wth)
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
    if i==0:
        ax.set_ylabel('Intra-chain Contact',fontsize=size)
    # ax.set_xscale('log')
    ax.set_xlim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

    ax=axs[1,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.scatter(zabs_tp_mean[j,i],p_inter[j,i],color=rgb,s=20*wth)
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
    ax.set_xlabel('Radial Position ($\\mathrm{\\AA}$)\nat %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Inter-chain Contact',fontsize=size)
    # ax.set_xscale('log')
    ax.set_xlim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

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

plt.savefig(f'{figname}_cont_cor_1.png',format='png')
plt.savefig(f'{figname}_cont_cor_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nmd,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
typeidx=np.arange(ntp)

for i in range(nmd):
    ax=axs[0,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.scatter(zabs_tp_mean[i,j],p_intra[i,j],color=rgb,s=20*wth)
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
    if i==0:
        ax.set_ylabel('Intra-chain Contact',fontsize=size)
    # ax.set_xscale('log')
    ax.set_xlim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

    ax=axs[1,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.scatter(zabs_tp_mean[i,j],p_inter[i,j],color=rgb,s=20*wth)
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
    ax.set_xlabel('Radial Position ($\\mathrm{\\AA}$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Inter-chain Contact',fontsize=size)
    # ax.set_xscale('log')
    ax.set_xlim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nT)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Temperature (K)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nT),1-1/(2*nT),nT))
cbar.set_ticklabels([int(Ts[i]) for i in range(nT)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_cont_cor_2.png',format='png')
plt.savefig(f'{figname}_cont_cor_2.pdf',format='pdf')
# plt.show()
