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
dir='ori'
dir2='../dst/dst'
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
nhi_z=250
nhi_tht=50

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

z=np.zeros((nmd,nT,nhi_z))
tht=np.zeros((nmd,nT,nhi_tht))
p=np.zeros((nmd,nT,nhi_z,nhi_tht))
tht_mean=np.zeros((nmd,nT,nhi_z))
tht_std=np.zeros((nmd,nT,nhi_z))
thtabs=np.zeros((nmd,nT,nhi_tht))
pabs=np.zeros((nmd,nT,nhi_z,nhi_tht))
thtabs_mean=np.zeros((nmd,nT,nhi_z))
thtabs_std=np.zeros((nmd,nT,nhi_z))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        tht[i,j]=q['tht']
        p[i,j]=q['p']
        tht_mean[i,j]=q['tht_mean']
        tht_std[i,j]=q['tht_std']
        thtabs[i,j]=q['thtabs']
        pabs[i,j]=q['pabs']
        thtabs_mean[i,j]=q['thtabs_mean']
        thtabs_std[i,j]=q['thtabs_std']
tht*=180/np.pi
p/=180/np.pi
tht_mean*=180/np.pi
tht_std*=180/np.pi
thtabs*=180/np.pi
pabs/=180/np.pi
thtabs_mean*=180/np.pi
thtabs_std*=180/np.pi
tht_mean[tht_mean==0.0]=np.nan
tht_std[tht_std==0.0]=np.nan
thtabs_mean[tht_mean==0.0]=np.nan
thtabs_std[tht_std==0.0]=np.nan

pedge=np.zeros((nmd,nT,3,nhi_tht))
pabsedge=np.zeros((nmd,nT,3,nhi_tht))
pedge[:,:,0]=np.sum(p[:,:,round(nhi_z/2)-1:round(nhi_z/2)+1],axis=2)
pabsedge[:,:,0]=np.sum(pabs[:,:,round(nhi_z/2)-1:round(nhi_z/2)+1],axis=2)
z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir2}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]
        for k in range(nhi_z):
            if z1[i,j]+d[i,j]<z[i,j,k]<z2[i,j]-d[i,j]:
                pedge[i,j,1]+=p[i,j,k]
                pabsedge[i,j,1]+=pabs[i,j,k]
            if z1[i,j]-d[i,j]<z[i,j,k]<z1[i,j]+d[i,j] or z2[i,j]-d[i,j]<z[i,j,k]<z2[i,j]+d[i,j]:
                pedge[i,j,2]+=p[i,j,k]
                pabsedge[i,j,2]+=pabs[i,j,k]
pedge/=np.sum(pedge*180/nhi_tht,axis=-1,keepdims=True)
pabsedge/=np.sum(pabsedge*90/nhi_tht,axis=-1,keepdims=True)

pedge_rand=np.pi/360*np.sin(np.pi*tht[0,0]/180)
dpedge=2*np.arccos(np.clip(np.sum(np.sqrt(pedge)*np.sqrt(pedge_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_tht,axis=-1),-1.0,1.0))
# dpedge=np.sum((np.sqrt(pedge)-np.sqrt(pedge_rand[np.newaxis,np.newaxis,np.newaxis]))**2*180/nhi_tht,axis=-1)
pabsedge_rand=2*np.pi/360*np.sin(np.pi*thtabs[0,0]/180)
dpabsedge=2*np.arccos(np.clip(np.sum(np.sqrt(pabsedge)*np.sqrt(pabsedge_rand[np.newaxis,np.newaxis,np.newaxis])*90/nhi_tht,axis=-1),-1.0,1.0))
# dpabsedge=np.sum((np.sqrt(pabsedge)-np.sqrt(pabsedge_rand[np.newaxis,np.newaxis,np.newaxis]))**2*90/nhi_tht,axis=-1)

pslc=p/np.sum(p*180/nhi_tht,axis=-1,keepdims=True)
pabsslc=p/np.sum(pabs*90/nhi_tht,axis=-1,keepdims=True)
dp=2*np.arccos(np.clip(np.sum(np.sqrt(pslc)*np.sqrt(pedge_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_tht,axis=-1),-1.0,1.0))
dpabs=2*np.arccos(np.clip(np.sum(np.sqrt(pabsslc)*np.sqrt(pabsedge_rand[np.newaxis,np.newaxis,np.newaxis])*90/nhi_tht,axis=-1),-1.0,1.0))

"""Plot"""
# fig,axs=plt.subplots(nmd,nT,figsize=(15,12))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# # nodes=[0.00,1/4,2/4,3/4,1.00]
# color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# # cmap=plt.cm.rainbow
# norm=colors.LogNorm(vmin=1e-5,vmax=np.max(p))

# for i in range(nmd):
#     for j in range(nT):
#         ax=axs[i,j]
#         ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),tht[i,j],p[i,j].T,cmap=cmap,norm=norm)
#         ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),tht_mean[i,j],'c',linewidth=wth)
#         ax.autoscale()
#         ax.minorticks_on()
#         ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#         ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.yaxis.set_major_locator(MultipleLocator(45))
#         ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.spines['bottom'].set_linewidth(wth)
#         ax.spines['top'].set_linewidth(wth)
#         ax.spines['left'].set_linewidth(wth)
#         ax.spines['right'].set_linewidth(wth)
#         ax.set_xlim(-1,1)
#         if i==nmd-1:
#             ax.set_xlabel('Scaled Z\nat %d K'%Ts[j],fontsize=size)
#         if j==0:
#             ax.set_ylabel('Angle ($\\circ$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Probability Density',fontsize=size)
# # cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
# cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}.png',format='png')
# plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

# fig,axs=plt.subplots(3,nT,figsize=(15,12))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# ylabels=['Center','Inside','Surface']

# for i in range(nT):
#     for j in range(3):
#         ax=axs[j,i]
#         for k in range(nmd):
#             rgb=cmap(k/(nmd-1))
#             ax.plot(tht[k,i],pedge[k,i,j]*100,color=rgb,linewidth=wth)
#         ax.plot(tht[0,0],pedge_rand*100,'k--',linewidth=wth,label='Random')
#         ax.autoscale()
#         ax.minorticks_on()
#         ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#         ax.xaxis.set_major_locator(MultipleLocator(45))
#         ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.spines['bottom'].set_linewidth(wth)
#         ax.spines['top'].set_linewidth(wth)
#         ax.spines['left'].set_linewidth(wth)
#         ax.spines['right'].set_linewidth(wth)
#         if j==2:
#             ax.set_xlabel('Angle ($\\circ$)\nat %i K'%Ts[i],fontsize=size)
#         if i==0:
#             ax.set_ylabel('Probability Density\n($\\times10^{-2}$) at %s'%ylabels[j],fontsize=size)
#         if (i,j)==(0,0):
#             ax.legend(loc='best',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Structure Fluctuation',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
# cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_1.png',format='png')
# plt.savefig(f'{figname}_1.pdf',format='pdf')
# plt.show()

# fig,axs=plt.subplots(3,nmd,figsize=(15,12))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# color=['#087F5B','#0CA678','#20C997','#63E6BE']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# ylabels=['Center','Inside','Surface']

# for i in range(nmd):
#     for j in range(3):
#         ax=axs[j,i]
#         for k in range(nT):
#             rgb=cmap(k/(nT-1))
#             ax.plot(tht[i,k],pedge[i,k,j]*100,color=rgb,linewidth=wth)
#         ax.plot(tht[0,0],pedge_rand*100,'k--',linewidth=wth,label='Random')
#         ax.autoscale()
#         ax.minorticks_on()
#         ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#         ax.xaxis.set_major_locator(MultipleLocator(45))
#         ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#         ax.spines['bottom'].set_linewidth(wth)
#         ax.spines['top'].set_linewidth(wth)
#         ax.spines['left'].set_linewidth(wth)
#         ax.spines['right'].set_linewidth(wth)
#         if j==2:
#             ax.set_xlabel('Angle ($\\circ$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#         if i==0:
#             ax.set_ylabel('Probability Density\n($\\times10^{-2}$) at %s'%ylabels[j],fontsize=size)
#         if (i,j)==(0,0):
#             ax.legend(loc='best',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nT)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Temperature (K)',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nT),1-1/(2*nT),nT))
# cbar.set_ticklabels([int(Ts[i]) for i in range(nT)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_2.png',format='png')
# plt.savefig(f'{figname}_2.pdf',format='pdf')
# plt.show()

# fig,axs=plt.subplots(1,3,figsize=(15,5))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# xlabels=['Center','Inside','Surface']

# spc=3
# for i in range(3):
#     ax=axs[i]
#     for j in range(nmd):
#         rgb=cmap(j/(nmd-1))
#         ax.bar(np.arange(nT)*(nmd+1+spc)+j,dpedge[j,:,i],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
#     ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xticks(np.arange(nT)*(nmd+1+spc)+nmd/2,np.array(Ts).astype(int))
#     ax.set_xlabel(xlabels[i],fontsize=size)
#     ax.set_ylabel('Distance of Probability',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[-1].get_position().y1
# bot=axs[-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Structure Fluctuation',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
# cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_dist_1.png',format='png')
# plt.savefig(f'{figname}_dist_1.pdf',format='pdf')
# plt.show()

# fig,axs=plt.subplots(1,3,figsize=(15,5))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# color=['#087F5B','#0CA678','#20C997','#63E6BE']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# xlabels=['Center','Inside','Surface']

# spc=3
# for i in range(3):
#     ax=axs[i]
#     for j in range(nT):
#         rgb=cmap(j/(nT-1))
#         ax.bar(np.arange(nmd)*(nT+1+spc)+j,dpedge[:,j,i],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
#     ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xticks(np.arange(nmd)*(nT+1+spc)+nT/2,np.around(10**rg_flu,2))
#     ax.set_xlabel(xlabels[i],fontsize=size)
#     ax.set_ylabel('Distance of Probability',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[-1].get_position().y1
# bot=axs[-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nT)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Temperature (K)',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nT),1-1/(2*nT),nT))
# cbar.set_ticklabels([int(Ts[i]) for i in range(nT)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_dist_2.png',format='png')
# plt.savefig(f'{figname}_dist_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nT,figsize=(15,5))
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

for i in range(nT):
    ax=axs[i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),dp[j,i],color=rgb,linewidth=wth)
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
    ax.set_xlim(-1,1)
    ax.set_xlabel('Scaled Z at %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Distance of Probability',fontsize=size)

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

plt.savefig(f'{figname}_dist_1.png',format='png')
plt.savefig(f'{figname}_dist_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(15,5))
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
xlabels=['Center','Inside','Surface']

for i in range(nmd):
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),dp[i,j],color=rgb,linewidth=wth)
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
    ax.set_xlim(-1,1)
    ax.set_xlabel('Scaled Z\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Distance of Probability',fontsize=size)

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

plt.savefig(f'{figname}_dist_2.png',format='png')
plt.savefig(f'{figname}_dist_2.pdf',format='pdf')
plt.show()

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
norm=colors.LogNorm(vmin=1e-5,vmax=np.max(p))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),thtabs[i,j],pabs[i,j].T,cmap=cmap,norm=norm)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),thtabs_mean[i,j],'c',linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(45))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('Scaled Z\nat %d K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Angle ($\\circ$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Probability Density',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_abs.png',format='png')
plt.savefig(f'{figname}_abs.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(3,nT,figsize=(15,12))
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
ylabels=['Center','Inside','Surface']

for i in range(nT):
    for j in range(3):
        ax=axs[j,i]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(thtabs[k,i],pabsedge[k,i,j]*10,color=rgb,linewidth=wth)
        ax.plot(thtabs[0,0],pabsedge_rand*10,'k--',linewidth=wth,label='Random')
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.xaxis.set_major_locator(MultipleLocator(45))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        if j==2:
            ax.set_xlabel('Angle ($\\circ$)\nat %i K'%Ts[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Probability Density\n($\\times10^{-1}$) at %s'%ylabels[j],fontsize=size)
        if (i,j)==(0,0):
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
cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_abs_1.png',format='png')
plt.savefig(f'{figname}_abs_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(3,nmd,figsize=(15,12))
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
ylabels=['Center','Inside','Surface']

for i in range(nmd):
    for j in range(3):
        ax=axs[j,i]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(thtabs[i,k],pabsedge[i,k,j]*10,color=rgb,linewidth=wth)
        ax.plot(thtabs[0,0],pabsedge_rand*10,'k--',linewidth=wth,label='Random')
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.xaxis.set_major_locator(MultipleLocator(45))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        if j==2:
            ax.set_xlabel('Angle ($\\circ$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Probability Density\n($\\times10^{-1}$) at %s'%ylabels[j],fontsize=size)
        if (i,j)==(0,0):
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

plt.savefig(f'{figname}_abs_2.png',format='png')
plt.savefig(f'{figname}_abs_2.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,3,figsize=(15,5))
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
xlabels=['Center','Inside','Surface']

spc=3
for i in range(3):
    ax=axs[i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.bar(np.arange(nT)*(nmd+1+spc)+j,dpabsedge[j,:,i],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
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
    ax.set_xticks(np.arange(nT)*(nmd+1+spc)+nmd/2,np.array(Ts).astype(int))
    ax.set_xlabel(xlabels[i],fontsize=size)
    ax.set_ylabel('Distance of Probability',fontsize=size)

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

plt.savefig(f'{figname}_abs_dist_1.png',format='png')
plt.savefig(f'{figname}_abs_dist_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,3,figsize=(15,5))
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
xlabels=['Center','Inside','Surface']

spc=3
for i in range(3):
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.bar(np.arange(nmd)*(nT+1+spc)+j,dpabsedge[:,j,i],width=1.0,color=rgb,edgecolor='k',linewidth=wth)
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
    ax.set_xticks(np.arange(nmd)*(nT+1+spc)+nT/2,np.around(10**rg_flu,2))
    ax.set_xlabel(xlabels[i],fontsize=size)
    ax.set_ylabel('Distance of Probability',fontsize=size)

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

plt.savefig(f'{figname}_abs_dist_2.png',format='png')
plt.savefig(f'{figname}_abs_dist_2.pdf',format='pdf')
plt.show()