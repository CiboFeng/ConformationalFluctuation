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
dir='mag'
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
nhi_mag=50
maglim=[[-1,1],[-1,1],[-1,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
nq=len(maglim)

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
mag=np.zeros((nmd,nT,nhi_mag,nq))
p=np.zeros((nmd,nT,nhi_z,nhi_mag,nq))
mag_mean=np.zeros((nmd,nT,nhi_z,nq))
mag_std=np.zeros((nmd,nT,nhi_z,nq))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        mag[i,j]=q['Mag']
        p[i,j]=q['p']
        mag_mean[i,j]=q['mag_mean']
        mag_std[i,j]=q['mag_std']
mag_mean[mag_mean==0.0]=np.nan
mag_std[mag_std==0.0]=np.nan

z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir2}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]

"""Plot"""
figlabels=['x','y','z','abs_mean_x','abs_mean_y','abs_mean_z','abs_x','abs_y','abs_z','all']
ylabels=['X','Y','Z','abs(mean(X))','abs(mean(Y))','abs(mean(Z))','abs(X)','abs(Y)','abs(Z)','All']
# for i in range(nq):
#     fig,axs=plt.subplots(nmd,nT,figsize=(15,12))
#     wth=2
#     size=20
#     lenmaj=15
#     lenmin=8
#     xtick=0.1
#     ytick=20
#     figwth=0.87
#     cbarwth=0.03
#     # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
#     # nodes=[0.00,1/4,2/4,3/4,1.00]
#     color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
#     nodes=[0/3,1/3,2/3,3/3]
#     cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
#     # cmap=plt.cm.rainbow
#     norm=colors.Normalize(vmin=0,vmax=1)
#     # norm=colors.Normalize(vmin=0,vmax=1)

#     for j in range(nmd):
#         for k in range(nT):
#             ax=axs[j,k]
#             ax.pcolormesh(z[j,k]/((z2[j,k]-z1[j,k])/2+d[j,k]),mag[j,k,:,i],p[j,k,:,:,i].T,cmap=cmap,norm=norm)
#             ax.plot(z[j,k]/((z2[j,k]-z1[j,k])/2+d[j,k]),mag_mean[j,k,:,i],'c',linewidth=wth)
#             ax.autoscale()
#             ax.minorticks_on()
#             ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#             ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#             ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#             # ax.yaxis.set_major_locator(MultipleLocator(45))
#             ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#             ax.spines['bottom'].set_linewidth(wth)
#             ax.spines['top'].set_linewidth(wth)
#             ax.spines['left'].set_linewidth(wth)
#             ax.spines['right'].set_linewidth(wth)
#             ax.set_xlim(-1,1)
#             if j==nmd-1:
#                 ax.set_xlabel('Scaled Z\nat %d K'%Ts[k],fontsize=size)
#             if k==0:
#                 ax.set_ylabel('Magnetization of %s\nof Flu.=%.2f'%(ylabels[i],10**rg_flu[j]),fontsize=size)

#     fig.subplots_adjust(right=figwth)
#     plt.tight_layout(rect=[0,0,figwth,1])
#     top=axs[0,-1].get_position().y1
#     bot=axs[-1,-1].get_position().y0
#     cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
#     sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
#     sm.set_array([])
#     cbar=plt.colorbar(sm,cax=cbar_ax)
#     cbar.set_label('Probability Density',fontsize=size)
#     # cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
#     # cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
#     cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
#     cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
#     cbar.outline.set_linewidth(wth)

#     plt.savefig(f'{figname}_{figlabels[i]}.png',format='png')
#     plt.savefig(f'{figname}_{figlabels[i]}.pdf',format='pdf')
#     plt.show()

fig,axs=plt.subplots(4,nT,figsize=(15,12))
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
ylabels=['X','Y','Z','All']

for i in range(nT):
    for j in range(4):
        ax=axs[j,i]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(z[k,i]/((z2[k,i]-z1[k,i])/2+d[k,i]),mag_mean[k,i,:,j+3],color=rgb,linestyle='--',linewidth=wth)
            ax.plot(z[k,i]/((z2[k,i]-z1[k,i])/2+d[k,i]),mag_mean[k,i,:,j+6],color=rgb,linewidth=wth)
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
        ax.set_xlim(-1,1)
        ax.set_ylim(0,1)
        if j==3:
            ax.set_xlabel('Scaled Z\nat %i K'%Ts[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Magnetization\nof %s'%ylabels[j],fontsize=size)

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

plt.savefig(f'{figname}_1.png',format='png')
plt.savefig(f'{figname}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(4,nmd,figsize=(15,12))
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
ylabels=['X','Y','Z','All']

for i in range(nmd):
    for j in range(4):
        ax=axs[j,i]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(z[i,k]/((z2[i,k]-z1[i,k])/2+d[i,k]),mag_mean[i,k,:,j+3],color=rgb,linestyle='--',linewidth=wth)
            ax.plot(z[i,k]/((z2[i,k]-z1[i,k])/2+d[i,k]),mag_mean[i,k,:,j+6],color=rgb,linewidth=wth)
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
        ax.set_xlim(-1,1)
        ax.set_ylim(0,1)
        if j==3:
            ax.set_xlabel('Scaled Z\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        if i==0:
            ax.set_ylabel('Magnetization\nof %s'%ylabels[j],fontsize=size)

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