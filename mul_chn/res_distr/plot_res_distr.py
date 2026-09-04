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
from pyw import pyw

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir1='../dst/dst'
dir='res_distr'
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

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(15,12))
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
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=np.max(pz_seq))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        # ax.imshow(p_seq[i,j],cmap=cmap)
        ax.pcolormesh(z_seq[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),np.arange(nat),pz_seq[i,j],cmap=cmap)
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
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('Scaled Z\nat %d K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Residue Index\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Ratio of Each Residue',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_seq.png',format='png')
plt.savefig(f'{figname}_seq.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nT,figsize=(15,4))
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
seqidx=np.arange(nat)

for i in range(nT):
    ax=axs[i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(zabs_seq_mean[j,i],color=rgb,linewidth=wth)
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
    ax.set_xlabel('Residue Index\nat %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radial Position ($\\mathrm{\\AA}$)',fontsize=size)
    # ax.set_xscale('log')
    ax.set_ylim(1.05*np.min(zabs_seq_mean)-0.05*np.max(zabs_seq_mean),-0.05*np.min(zabs_seq_mean)+1.05*np.max(zabs_seq_mean))

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

plt.savefig(f'{figname}_seq_mean_1.png',format='png')
plt.savefig(f'{figname}_seq_mean_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(15,4))
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
seqidx=np.arange(nat)

for i in range(nmd):
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(zabs_seq_mean[i,j],color=rgb,linewidth=wth)
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
    ax.set_xlabel('Residue Index\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radial Position ($\\mathrm{\\AA}$)',fontsize=size)
    # ax.set_xscale('log')
    ax.set_ylim(1.05*np.min(zabs_seq_mean)-0.05*np.max(zabs_seq_mean),-0.05*np.min(zabs_seq_mean)+1.05*np.max(zabs_seq_mean))

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
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

plt.savefig(f'{figname}_seq_mean_2.png',format='png')
plt.savefig(f'{figname}_seq_mean_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nT,figsize=(15,12))
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
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=np.max(pz_tp))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        # ax.imshow(p_tp[i,j],cmap=cmap)
        ax.pcolormesh(z_tp[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),np.arange(ntp),pz_tp[i,j],cmap=cmap)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(-1,1)
        ax.set_yticks(np.arange(ntp),res_sort)
        if i==nmd-1:
            ax.set_xlabel('Scaled Z\nat %d K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Ratio of Each Residue Type',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_type.png',format='png')
plt.savefig(f'{figname}_type.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nT,figsize=(15,4))
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
    ax=axs[i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.scatter(np.arange(ntp),zabs_tp_mean[j,i],color=rgb,s=20*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_xlabel('%i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radial Position ($\\mathrm{\\AA}$)',fontsize=size)
    # ax.set_xscale('log')
    ax.set_ylim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

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

plt.savefig(f'{figname}_type_mean_1.png',format='png')
plt.savefig(f'{figname}_type_mean_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(15,4))
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
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.scatter(np.arange(ntp),zabs_tp_mean[i,j],color=rgb,s=20*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_xlabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radial Position ($\\mathrm{\\AA}$)',fontsize=size)
    # ax.set_xscale('log')
    ax.set_ylim(1.05*np.nanmin(zabs_tp_mean)-0.05*np.nanmax(zabs_tp_mean),-0.05*np.nanmin(zabs_tp_mean)+1.05*np.nanmax(zabs_tp_mean))

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
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

plt.savefig(f'{figname}_type_mean_2.png',format='png')
plt.savefig(f'{figname}_type_mean_2.pdf',format='pdf')
# plt.show()
