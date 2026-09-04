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
dir1='../dst/dst'
dir='rg_slab'
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
nhirg0=200
nhiz=250
nhirg=50
rglim=[10.0,60.0]
nsms=10

"""Read Data"""
rg0=np.zeros((nmd,nT,nhirg0))
p_rg0=np.zeros((nmd,nT,nhirg0))
rg0_mean=np.zeros((nmd,nT))
rg0_std=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Probability')
        rg0[i,j]=q[0]
        p_rg0[i,j]=q[1]
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Mean')
        rg0_mean[i,j]=q[0,0]
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Deviation')
        rg0_std[i,j]=q[0,0]
rg0_flu=np.round(np.log10(rg0_std[:,0]/rg0_mean[:,0]),1)

z=np.zeros((nmd,nT,nhiz))
rg=np.zeros((nmd,nT,nhirg))
p=np.zeros((nmd,nT,nhiz,nhirg))
rg_mean=np.zeros((nmd,nT,nhiz))
rg_std=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        rg[i,j]=q['rg']
        p[i,j]=q['p']
        rg_mean[i,j]=q['rg_mean']
        rg_std[i,j]=q['rg_std']

pz=np.sum(p,axis=-1)
pz/=np.sum(pz*(np.max(z)-np.min(z))/nhiz,axis=-1,keepdims=True)
z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
pedge=np.zeros((nmd,nT,3,nhirg))
rgedge_mean=np.zeros((nmd,nT,3))
rgedge_std=np.zeros((nmd,nT,3))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]
        idxz=[[],[],[]]
        idxz[0]+=[round(nhiz/2)-1,round(nhiz/2)+1]
        for k in range(nhiz):
            if z1[i,j]+d[i,j]<z[i,j,k]<z2[i,j]-d[i,j]:
                idxz[1]+=[k]
            if z1[i,j]-d[i,j]<z[i,j,k]<z1[i,j]+d[i,j] or z2[i,j]-d[i,j]<z[i,j,k]<z2[i,j]+d[i,j]:
                idxz[2]+=[k]
        for k in range(len(idxz)):
            idx=(i,j,idxz[k])
            pedge[i,j,k]=np.sum(p[idx],axis=0)
            rgedge_mean[i,j,k]=np.sum(pz[idx]*rg_mean[idx])/np.sum(pz[idx])
            rgedge_std[i,j,k]=np.sqrt(np.sum(pz[idx]*rg_std[idx]**2)/np.sum(pz[idx])+ \
                    np.sum(pz[idx]*(rg_mean[idx]-rgedge_mean[i,j,k])**2)/np.sum(pz[idx]))
pedge/=np.sum(pedge*(np.max(rg)-np.min(rg))/nhirg,axis=-1,keepdims=True)

rg_mean[rg_mean==0.0]=np.nan
rg_std[rg_std==0.0]=np.nan

rg_mean_sms=np.zeros((nmd,nT,nhiz))
rg_std_sms=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        for k in range(nhiz):
            idx=(i,j,list(np.arange(max(0,k-round(nsms/2)),min(nhiz,k+round(nsms/2)))))
            rg_mean_sms[i,j,k]=np.nansum(pz[idx]*rg_mean[idx])/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0))
            rg_std_sms[i,j,k]=np.sqrt(np.nansum(pz[idx]*rg_std[idx]**2)/np.nansum(pz[idx]*np.heaviside(rg_std[idx],0.0))+ \
                    np.nansum(pz[idx]*(rg_mean[idx]-rg_mean_sms[i,j,k])**2)/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0)))

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
norm=colors.LogNorm(vmin=1e-4,vmax=np.max(p))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg[i,j],p[i,j].T,cmap=cmap)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_mean[i,j],'c',linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(45))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('Scaled Z\nat %d K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)\nof Flu.=%.2f'%10**rg0_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Probability Density',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(4,nT+1,figsize=(20,12))
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
ylabels=['Center','Inside','Surface','Isolated']

for i in range(nT):
    for j in range(3):
        ax=axs[j,i]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(rg[k,i],pedge[k,i,j],color=rgb,linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax.xaxis.set_major_locator(MultipleLocator(45))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(1.05*rglim[0]-0.05*rglim[1],-0.05*rglim[0]+1.05*rglim[1])
        # if j==2:
        #     ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)\nof Flu.=%.2f'%10**rg0_flu[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Probability Density\nat %s'%ylabels[j],fontsize=size)

for i in range(nT):
    ax=axs[3,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(rg0[j,i],p_rg0[j,i],color=rgb,linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(45))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(1.05*rglim[0]-0.05*rglim[1],-0.05*rglim[0]+1.05*rglim[1])
    ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)\nat %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Probability Density\nat %s'%ylabels[j],fontsize=size)

spc=3
for i in range(3):
    ax=axs[i,nT]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.errorbar(np.arange(nT)*(nmd+1+spc)+j,rgedge_mean[j,:,i],yerr=rgedge_std[j,:,i],
                        linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                        ecolor=rgb,elinewidth=wth,capsize=2*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_ylim(10,42)
    ax.set_xticks(np.arange(nT)*(nmd+1+spc)+nmd/2,np.array(Ts).astype(int))
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)

ax=axs[3,nT]
for j in range(nmd):
    rgb=cmap(j/(nmd-1))
    ax.errorbar(np.arange(nT)*(nmd+1+spc)+j,rg0_mean[j,:],yerr=rg0_std[j,:],
                    linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                    ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_ylim(10,42)
ax.set_xticks(np.arange(nT)*(nmd+1+spc)+nmd/2,np.array(Ts).astype(int))
ax.set_xlabel('Temperature',fontsize=size)
ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
if i==0:
    ax.legend(loc='best',fontsize=size)

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

plt.savefig(f'{figname}_1.png',format='png')
plt.savefig(f'{figname}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(4,nmd+1,figsize=(15,12))
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
ylabels=['Center','Inside','Surface','Isolated']

for i in range(nmd):
    for j in range(3):
        ax=axs[j,i]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(rg[i,k],pedge[i,k,j],color=rgb,linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax.xaxis.set_major_locator(MultipleLocator(45))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(1.05*rglim[0]-0.05*rglim[1],-0.05*rglim[0]+1.05*rglim[1])
        # if j==2:
        #     ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)\nof Flu.=%.2f'%10**rg0_flu[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Probability Density\nat %s'%ylabels[j],fontsize=size)

for i in range(nmd):
    ax=axs[3,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(rg0[i,j],p_rg0[i,j],color=rgb,linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(45))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(1.05*rglim[0]-0.05*rglim[1],-0.05*rglim[0]+1.05*rglim[1])
    ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)\nof Flu.=%.2f'%10**rg0_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Probability Density\nat %s'%ylabels[j],fontsize=size)

spc=3
for i in range(3):
    ax=axs[i,nmd]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.errorbar(np.arange(nmd)*(nT+1+spc)+j,rgedge_mean[:,j,i],yerr=rgedge_std[:,j,i],
                        linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                        ecolor=rgb,elinewidth=wth,capsize=2*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_ylim(10,42)
    ax.set_xticks(np.arange(nmd)*(nT+1+spc)+nT/2,np.around(10**rg0_flu,2))
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)

ax=axs[3,nmd]
for j in range(nT):
    rgb=cmap(j/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+1+spc)+j,rg0_mean[:,j],yerr=rg0_std[:,j],
                    linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                    ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_ylim(10,42)
ax.set_xticks(np.arange(nmd)*(nT+1+spc)+nT/2,np.around(10**rg0_flu,2))
ax.set_xlabel('Structure Fluctuation',fontsize=size)
ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
if i==0:
    ax.legend(loc='best',fontsize=size)

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

plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()

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
xlabels=['Mean','Dev.','Flu.']

for i in range(nT):
    ax=axs[0,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        # ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),rg_std_sms[j,i]/rg_mean_sms[j,i],color=rgb,linewidth=wth)
        ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),rg_std[j,i]/rg_mean[j,i],color=rgb,linewidth=wth)
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
    # ax.set_xlim(-275,275)
    ax.set_xlim(-1,1)
    ax.set_ylim(0.195,0.245)
    ax.set_xlabel('Scaled Z',fontsize=size)
    if i==0:
        ax.set_ylabel('Radius of Gyration\nFluctuation',fontsize=size)

    ax=axs[1,i]
    spc=3
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        y=rgedge_mean[j,i,1]/np.mean(rgedge_mean[:,:,1])-1
        ax.bar([0*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
        y=rgedge_std[j,i,1]/np.mean(rgedge_std[:,:,1])-1
        ax.bar([1*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
        y=rgedge_std[j,i,1]/rgedge_mean[j,i,1]/np.mean(rgedge_std[:,:,1]/rgedge_mean[:,:,1])-1
        ax.bar([2*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_ylim(-0.18,0.15)
    # ax.set_yscale('log')
    ax.set_xticks(np.arange(3)*(nmd+spc)+nmd/2,xlabels)
    ax.set_xlabel('%i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radius of Gyration\nFluctuation',fontsize=size)

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

plt.savefig(f'{figname}_flu_1.png',format='png')
plt.savefig(f'{figname}_flu_1.pdf',format='pdf')
plt.show()

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
xlabels=['Mean','Dev.','Flu.']

for i in range(nmd):
    ax=axs[0,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        # ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_std_sms[i,j]/rg_mean_sms[i,j],color=rgb,linewidth=wth)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_std[i,j]/rg_mean[i,j],color=rgb,linewidth=wth)
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
    # ax.set_xlim(-275,275)
    ax.set_xlim(-1,1)
    ax.set_ylim(0.195,0.245)
    ax.set_xlabel('Scaled Z',fontsize=size)
    if i==0:
        ax.set_ylabel('Radius of Gyration\nFluctuation',fontsize=size)

    ax=axs[1,i]
    spc=3
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        y=rgedge_mean[i,j,1]/np.mean(rgedge_mean[:,:,1])-1
        ax.bar([0*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
        y=rgedge_std[i,j,1]/np.mean(rgedge_std[:,:,1])-1
        ax.bar([1*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
        y=rgedge_std[i,j,1]/rgedge_mean[i,j,1]/np.mean(rgedge_std[:,:,1]/rgedge_mean[:,:,1])-1
        ax.bar([2*(nT+spc)+j+1],[y],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_ylim(-0.18,0.15)
    # ax.set_yscale('log')
    ax.set_xticks(np.arange(3)*(nT+spc)+nT/2,xlabels)
    ax.set_xlabel('Flu.=%.2f'%10**rg0_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Radius of Gyration\nFluctuation',fontsize=size)

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

plt.savefig(f'{figname}_flu_2.png',format='png')
plt.savefig(f'{figname}_flu_2.pdf',format='pdf')
plt.show()