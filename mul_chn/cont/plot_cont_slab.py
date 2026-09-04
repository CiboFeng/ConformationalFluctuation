"""Import Modules"""
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
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

""""""
dir0='../../max_entr/rg/rg_heat'
dir='cont_slab'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nfr=10000
nch=100
nrep=20
ntp=len(res)
dt=10e-6*10000 ### ns

""""""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

P0=np.zeros((nmd,nT,nat,nat))
Pintra=np.zeros((nmd,nT,nat,nat))
Pinter=np.zeros((nmd,nT,nat,nat))
Pchn=np.zeros((nmd,nT,nfr,nch,nch))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        Pintra[i,j]=q['Pintra']
        Pinter[i,j]=q['Pinter']
        Pchn[i,j]=q['Pchn']

Pintra_=Pintra.copy()
Pinter_=Pinter.copy()
for i in range(4):
    for j in range(nat-i):
        Pintra_[:,:,j,i+j]=0.0
        Pinter_[:,:,j,i+j]=0.0

Plintra=np.sum(Pintra_,axis=-1)
Plinter=np.sum(Pinter_,axis=-1)

Psintra=np.zeros((nmd,nT,nat))
Psinter=np.zeros((nmd,nT,nat))
for i in range(nat):
    for j in range(nat-i):
        Psintra[:,:,i]+=Pintra[:,:,i+j,j]
        Psinter[:,:,i]+=Pinter[:,:,i+j,j]
    Psintra[:,:,i]/=nat-i
    Psinter[:,:,i]/=nat-i

pintra=np.zeros((nmd,nT,ntp,ntp))
pinter=np.zeros((nmd,nT,ntp,ntp))
num=np.zeros((ntp,ntp),dtype=int)
for i in range(nat):
    for j in range(i+1):
        pintra[:,:,seq[i],seq[j]]+=Pintra_[:,:,i,j]
        pinter[:,:,seq[i],seq[j]]+=Pinter_[:,:,i,j]
        num[seq[i],seq[j]]+=1
        if seq[i]!=seq[j]:
            pintra[:,:,seq[j],seq[i]]+=Pintra_[:,:,i,j]
            pinter[:,:,seq[j],seq[i]]+=Pinter_[:,:,i,j]
            num[seq[j],seq[i]]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            pintra[:,:,i,j]/=num[i,j]
            pinter[:,:,i,j]/=num[i,j]
pintra[pintra==0.0]=np.nan
pinter[pinter==0.0]=np.nan

plintra=np.nansum(pintra,axis=-1)
plinter=np.nansum(pinter,axis=-1)
plintra[plintra==0.0]=np.nan
plinter[plinter==0.0]=np.nan

idx_sort=[res.index(i) for i in res_sort]
pintra=pintra[:,:,idx_sort]
pintra=pintra[:,:,:,idx_sort]
pinter=pinter[:,:,idx_sort]
pinter=pinter[:,:,:,idx_sort]
plintra=plintra[:,:,idx_sort]
plinter=plinter[:,:,idx_sort]

Psintra*=1e1
plintra*=1e1
Plinter*=1e1
plinter*=1e3

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.85
cbarwth=0.03
cbargap=0.05
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap0=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm0=colors.LogNorm(vmin=1e-3,vmax=1)
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-3,vmax=1)
# qshow=Pintra[-1:,:1]-Pintra
# xmin=-np.min(qshow)/np.max(np.abs(qshow))
# xmax=np.max(qshow)/np.max(np.abs(qshow))
# color=[(0,0,xmin),(1,1,1),(xmax,0,0)]
# nodes=[0,-np.min(qshow)/(np.max(qshow)-np.min(qshow)),1]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=np.min(qshow),vmax=np.max(qshow))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        if (i,j)==(nmd-1,0):
            ax.imshow(Pintra[i,j],cmap=cmap0,norm=norm0)
        else:
            ax.imshow(Pintra[-1,0]-Pintra[i,j],cmap=cmap,norm=norm)
        ax.invert_yaxis()
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
        # ax.set_xlim(-4,84)
        if i==nmd-1:
            ax.set_xlabel(f'{int(Ts[j])} K',fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbargap+cbarwth/2,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap0,norm=norm0)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_axes([figwth+cbargap,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('Contact Probability Difference',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar2.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar2.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar2.outline.set_linewidth(wth)
cbar2.ax.yaxis.set_ticks_position('left')
cbar2.ax.yaxis.set_label_position('left')
log=np.arange(np.ceil(np.log10(norm.vmin)),np.floor(np.log10(norm.vmax))+1)
cbar2.set_ticks([10**x for x in log])
cbar2.set_ticklabels(['$-10^{%s}$'%int(x) for x in log])
cbar2.ax.invert_yaxis()

plt.savefig(f'{figname}_seq_intra.png',format='png')
plt.savefig(f'{figname}_seq_intra.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nT,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.85
cbarwth=0.03
cbargap=0.05
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap0=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm0=colors.LogNorm(vmin=1e-4,vmax=1e-2)
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=1e-2)
# qshow=Pinter[-1:,:1]-Pinter
# xmin=-np.min(qshow)/np.max(np.abs(qshow))
# xmax=np.max(qshow)/np.max(np.abs(qshow))
# color=[(0,0,xmin),(1,1,1),(xmax,0,0)]
# nodes=[0,-np.min(qshow)/(np.max(qshow)-np.min(qshow)),1]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=np.min(qshow),vmax=np.max(qshow))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        if (i,j)==(nmd-1,0):
            ax.imshow(Pinter[i,j],cmap=cmap0,norm=norm0)
        else:
            ax.imshow(Pinter[-1,0]-Pinter[i,j],cmap=cmap,norm=norm)
        ax.invert_yaxis()
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
        # ax.set_xlim(-4,84)
        if i==nmd-1:
            ax.set_xlabel(f'{int(Ts[j])} K',fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbargap+cbarwth/2,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap0,norm=norm0)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_axes([figwth+cbargap,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('Contact Probability Difference',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar2.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar2.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar2.outline.set_linewidth(wth)
cbar2.ax.yaxis.set_ticks_position('left')
cbar2.ax.yaxis.set_label_position('left')
log=np.arange(np.ceil(np.log10(norm.vmin)),np.floor(np.log10(norm.vmax))+1)
cbar2.set_ticks([10**x for x in log])
cbar2.set_ticklabels(['$-10^{%s}$'%int(x) for x in log])
cbar2.ax.invert_yaxis()

plt.savefig(f'{figname}_seq_inter.png',format='png')
plt.savefig(f'{figname}_seq_inter.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nT,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.85
cbarwth=0.03
cbargap=0.05
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap0=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm0=colors.LogNorm(vmin=1e-4,vmax=1)
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=1)
# qshow=pintra[-1:,:1]-pintra
# xmin=-np.min(qshow)/np.max(np.abs(qshow))
# xmax=np.max(qshow)/np.max(np.abs(qshow))
# color=[(0,0,xmin),(1,1,1),(xmax,0,0)]
# nodes=[0,-np.min(qshow)/(np.max(qshow)-np.min(qshow)),1]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=np.min(qshow),vmax=np.max(qshow))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        if (i,j)==(nmd-1,0):
            ax.imshow(pintra[i,j],cmap=cmap0,norm=norm0)
        else:
            ax.imshow(pintra[-1,0]-pintra[i,j],cmap=cmap,norm=norm)
        ax.invert_yaxis()
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        # ax.set_xscale('log')
        # ax.set_xlim(-4,84)
        ax.set_xticks(np.arange(ntp),res_sort)
        ax.set_yticks(np.arange(ntp),res_sort)
        if i==nmd-1:
            ax.set_xlabel(f'{int(Ts[j])} K',fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbargap+cbarwth/2,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap0,norm=norm0)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_axes([figwth+cbargap,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('Contact Probability Difference',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar2.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar2.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar2.outline.set_linewidth(wth)
cbar2.ax.yaxis.set_ticks_position('left')
cbar2.ax.yaxis.set_label_position('left')
log=np.arange(np.ceil(np.log10(norm.vmin)),np.floor(np.log10(norm.vmax))+1)
cbar2.set_ticks([10**x for x in log])
cbar2.set_ticklabels(['$-10^{%s}$'%int(x) for x in log])
cbar2.ax.invert_yaxis()

plt.savefig(f'{figname}_type_intra.png',format='png')
plt.savefig(f'{figname}_type_intra.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nT,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.85
cbarwth=0.03
cbargap=0.05
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap0=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm0=colors.LogNorm(vmin=1e-4,vmax=1e-2)
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=1e-2)
# qshow=pinter[-1:,:1]-pinter
# xmin=-np.min(qshow)/np.max(np.abs(qshow))
# xmax=np.max(qshow)/np.max(np.abs(qshow))
# color=[(0,0,xmin),(1,1,1),(xmax,0,0)]
# nodes=[0,-np.min(qshow)/(np.max(qshow)-np.min(qshow)),1]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=np.min(qshow),vmax=np.max(qshow))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        if (i,j)==(nmd-1,0):
            ax.imshow(pinter[i,j],cmap=cmap0,norm=norm0)
        else:
            ax.imshow(pinter[-1,0]-pinter[i,j],cmap=cmap,norm=norm)
        ax.invert_yaxis()
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        # ax.set_xscale('log')
        # ax.set_xlim(-4,84)
        ax.set_xticks(np.arange(ntp),res_sort)
        ax.set_yticks(np.arange(ntp),res_sort)
        if i==nmd-1:
            ax.set_xlabel(f'{int(Ts[j])} K',fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbargap+cbarwth/2,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap0,norm=norm0)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_axes([figwth+cbargap,bot,cbarwth/2,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('Contact Probability Difference',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar2.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar2.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar2.outline.set_linewidth(wth)
cbar2.ax.yaxis.set_ticks_position('left')
cbar2.ax.yaxis.set_label_position('left')
log=np.arange(np.ceil(np.log10(norm.vmin)),np.floor(np.log10(norm.vmax))+1)
cbar2.set_ticks([10**x for x in log])
cbar2.set_ticklabels(['$-10^{%s}$'%int(x) for x in log])
cbar2.ax.invert_yaxis()

plt.savefig(f'{figname}_type_inter.png',format='png')
plt.savefig(f'{figname}_type_inter.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nT,3,figsize=(15,10))
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
# cmap=plt.cm.rainbow

for i in range(nT):
    ax=axs[i,0]
    if i==0:
        ax2=ax.twinx()
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        if (i,j)==(0,nmd-1):
            ax2.plot(np.arange(nat),Plintra[j,i],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Plintra[j,i],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Plintra[j,i]-Plintra[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Plintra-Plintra[-1:,:1])
    ymax=np.nanmax(Plintra-Plintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nT-1:
        ax.set_xlabel(f'Residue Index',fontsize=size)
    ax.set_ylabel('Difference\nat %i K'%Ts[i],fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum',fontsize=size)
    ax2.set_ylim(4,12)

    ax=axs[i,1]
    if i==0:
        ax2=ax.twinx()
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        if (i,j)==(0,nmd-1):
            ax2.plot(np.arange(nat),Psintra[j,i],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Psintra[j,i],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Psintra[j,i]-Psintra[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Psintra-Psintra[-1:,:1])
    ymax=np.nanmax(Psintra-Psintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nT-1:
        ax.set_xlabel(f'Sequence Distance',fontsize=size)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    # ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Mean ($\\times{10^{-1}}$)',fontsize=size)

    ax=axs[i,2]
    if i==0:
        ax2=ax.twinx()
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        if (i,j)==(0,nmd-1):
            ax2.scatter(np.arange(ntp),plintra[j,i],color='r',linewidth=4*wth,zorder=2)
            ax2.scatter(np.arange(ntp),plintra[j,i],color=rgb,linewidth=wth,zorder=2)
        else:
            ax.scatter(np.arange(ntp),plintra[j,i]-plintra[-1,0],color=rgb,linewidth=wth,zorder=2)
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis='x',which='major',linestyle='-',linewidth=wth,alpha=0.5,zorder=1)
    ax.grid(axis='y',visible=False)
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xscale('log')
    ymin=np.nanmin(plintra-plintra[-1:,:1])
    ymax=np.nanmax(plintra-plintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-1}}$)',fontsize=size)
    ax2.set_ylim(3,11)

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

plt.savefig(f'{figname}_sum_intra_1.png',format='png')
plt.savefig(f'{figname}_sum_intra_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,3,figsize=(15,10))
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
# cmap=plt.cm.rainbow

for i in range(nmd):
    ax=axs[i,0]
    if i==nmd-1:
        ax2=ax.twinx()
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        if (i,j)==(nmd-1,0):
            ax2.plot(np.arange(nat),Plintra[i,j],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Plintra[i,j],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Plintra[i,j]-Plintra[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Plintra-Plintra[-1:,:1])
    ymax=np.nanmax(Plintra-Plintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nmd-1:
        ax.set_xlabel(f'Residue Index',fontsize=size)
    ax.set_ylabel('Difference\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum',fontsize=size)
    ax2.set_ylim(4,12)

    ax=axs[i,1]
    if i==nmd-1:
        ax2=ax.twinx()
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        if (i,j)==(nmd-1,0):
            ax2.plot(np.arange(nat),Psintra[i,j],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Psintra[i,j],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Psintra[i,j]-Psintra[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Psintra-Psintra[-1:,:1])
    ymax=np.nanmax(Psintra-Psintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nmd-1:
        ax.set_xlabel(f'Sequence Distance',fontsize=size)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    # ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Mean ($\\times{10^{-1}}$)',fontsize=size)

    ax=axs[i,2]
    if i==nmd-1:
        ax2=ax.twinx()
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        if (i,j)==(nmd-1,0):
            ax2.scatter(np.arange(ntp),plintra[i,j],color='r',linewidth=4*wth,zorder=2)
            ax2.scatter(np.arange(ntp),plintra[i,j],color=rgb,linewidth=wth,zorder=2)
        else:
            ax.scatter(np.arange(ntp),plintra[i,j]-plintra[-1,0],color=rgb,linewidth=wth,zorder=2)
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis='x',which='major',linestyle='-',linewidth=wth,alpha=0.5,zorder=1)
    ax.grid(axis='y',visible=False)
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xscale('log')
    ymin=np.nanmin(plintra-plintra[-1:,:1])
    ymax=np.nanmax(plintra-plintra[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-1}}$)',fontsize=size)
    ax2.set_ylim(3,11)

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

plt.savefig(f'{figname}_sum_intra_2.png',format='png')
plt.savefig(f'{figname}_sum_intra_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nT,2,figsize=(15,10))
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
# cmap=plt.cm.rainbow

for i in range(nT):
    ax=axs[i,0]
    if i==0:
        ax2=ax.twinx()
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        if (i,j)==(0,nmd-1):
            ax2.plot(np.arange(nat),Plinter[j,i],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Plinter[j,i],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Plinter[j,i]-Plinter[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Plinter-Plinter[-1:,:1])
    ymax=np.nanmax(Plinter-Plinter[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nT-1:
        ax.set_xlabel(f'Residue Index',fontsize=size)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)\nat %i K'%Ts[i],fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-1}}$)',fontsize=size)
    ax2.set_ylim(0.5,3.5)

    ax=axs[i,1]
    if i==0:
        ax2=ax.twinx()
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        if (i,j)==(0,nmd-1):
            ax2.scatter(np.arange(ntp),plinter[j,i],color='r',linewidth=4*wth,zorder=2)
            ax2.scatter(np.arange(ntp),plinter[j,i],color=rgb,linewidth=wth,zorder=2)
        else:
            ax.scatter(np.arange(ntp),plinter[j,i]-plinter[-1,0],color=rgb,linewidth=wth,zorder=2)
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis='x',which='major',linestyle='-',linewidth=wth,alpha=0.5,zorder=1)
    ax.grid(axis='y',visible=False)
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xscale('log')
    ymin=np.nanmin(plinter-plinter[-1:,:1])
    ymax=np.nanmax(plinter-plinter[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_ylabel('Difference ($\\times{10^{-3}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_major_locator(MultipleLocator(4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-3}}$)',fontsize=size)
    ax2.set_ylim(3,14)

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

plt.savefig(f'{figname}_sum_inter_1.png',format='png')
plt.savefig(f'{figname}_sum_inter_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,2,figsize=(15,10))
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
# cmap=plt.cm.rainbow

for i in range(nmd):
    ax=axs[i,0]
    if i==nmd-1:
        ax2=ax.twinx()
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        if (i,j)==(nmd-1,0):
            ax2.plot(np.arange(nat),Plinter[i,j],color='r',linewidth=2*wth)
            ax2.plot(np.arange(nat),Plinter[i,j],color=rgb,linewidth=wth)
        else:
            ax.plot(np.arange(nat),Plinter[i,j]-Plinter[-1,0],color=rgb,linewidth=wth)
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
    ymin=np.nanmin(Plinter-Plinter[-1:,:1])
    ymax=np.nanmax(Plinter-Plinter[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    if i==nmd-1:
        ax.set_xlabel(f'Residue Index',fontsize=size)
    ax.set_ylabel('Difference ($\\times{10^{-1}}$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-1}}$)',fontsize=size)
    ax2.set_ylim(0.5,3.5)

    ax=axs[i,1]
    if i==nmd-1:
        ax2=ax.twinx()
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        if (i,j)==(nmd-1,0):
            ax2.scatter(np.arange(ntp),plinter[i,j],color='r',linewidth=4*wth,zorder=2)
            ax2.scatter(np.arange(ntp),plinter[i,j],color=rgb,linewidth=wth,zorder=2)
        else:
            ax.scatter(np.arange(ntp),plinter[i,j]-plinter[-1,0],color=rgb,linewidth=wth,zorder=2)
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmin,labelsize=0.8*size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis='x',which='major',linestyle='-',linewidth=wth,alpha=0.5,zorder=1)
    ax.grid(axis='y',visible=False)
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xscale('log')
    ymin=np.nanmin(plinter-plinter[-1:,:1])
    ymax=np.nanmax(plinter-plinter[-1:,:1])
    ax.set_ylim(1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin)
    ax.set_xticks(np.arange(ntp),res_sort)
    ax.set_ylabel('Difference ($\\times{10^{-3}}$)',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='r',labelcolor='r')
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='r',labelcolor='r')
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_major_locator(MultipleLocator(4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color('r')
    ax2.yaxis.label.set_color('r')
    ax2.set_ylabel('Contact Sum ($\\times{10^{-3}}$)',fontsize=size)
    ax2.set_ylim(3,14)

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

plt.savefig(f'{figname}_sum_inter_2.png',format='png')
plt.savefig(f'{figname}_sum_inter_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nT,figsize=(18,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-17,vmax=1e4)
# norm=colors.Normalize(vmin=0,vmax=1e-20)

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.imshow(Pchn[i,j,0],cmap=cmap,norm=norm)
        ax.invert_yaxis()
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
        # ax.set_xlim(-4,84)
        if i==nmd-1:
            ax.set_xlabel(f'{int(Ts[j])} K',fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_chn.png',format='png')
plt.savefig(f'{figname}_chn.pdf',format='pdf')
# plt.show()

