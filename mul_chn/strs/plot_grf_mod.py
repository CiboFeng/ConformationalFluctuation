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
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir='grf_mod'
figname1='relx_time'
figname2=dir
figname3='0_shear_mod'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
nomg=100

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

tau=np.zeros((nmd,nT,nfr,nch-1))
for i in range(nmd):
    for j in range(nT):
        tau[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npy')
tau_mean=np.mean(tau,axis=2)
# omg=np.logspace(3,7,nomg)
omg=np.linspace(0,2e6,nomg)
# u=omg[np.newaxis,np.newaxis,np.newaxis,np.newaxis]*tau[:,:,:,:,np.newaxis]
# G1=u**2/(1+u**2)
# G2=u/(1+u**2)+omg
# G1=np.mean(np.mean(G1,axis=-2),axis=-2)
# G2=np.mean(np.mean(G2,axis=-2),axis=-2)
# G1=np.zeros((nmd,nT,nomg))
# G2=np.zeros((nmd,nT,nomg))
# for i in range(nmd):
#     for j in range(nT):
#         u=omg[np.newaxis,np.newaxis]*tau[i,j,:,:,np.newaxis]
#         G1_tmp=u**2/(1+u**2)
#         G2_tmp=u/(1+u**2)
#         G1[i,j]=np.mean(np.mean(G1_tmp,axis=-2),axis=-2)
#         G2[i,j]=np.mean(np.mean(G2_tmp,axis=-2),axis=-2)
G1=np.zeros((nmd,nT,nomg))
G2=np.zeros((nmd,nT,nomg))
for i in range(nmd):
    for j in range(nT):
        u=omg[np.newaxis]*tau_mean[i,j,:,np.newaxis]
        G1_tmp=u**2/(1+u**2)
        G2_tmp=u/(1+u**2)
        G1[i,j]=np.mean(G1_tmp,axis=-2)
        G2[i,j]=np.mean(G2_tmp,axis=-2)

X01=np.mean(tau**2,axis=-1)/np.mean(tau,axis=-1)**2/np.array(Ts)[np.newaxis,:,np.newaxis]
X01_mean=np.mean(X01,axis=-1)
X01_std=np.std(X01,axis=-1)
X02=np.mean(tau,axis=-1)*np.array(Ts)[np.newaxis,:,np.newaxis]
X02_mean=np.mean(X02,axis=-1)
X02_std=np.std(X02,axis=-1)

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.88
cbarwth=0.03
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=100)
# yticks=[[[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4]],
#         [[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4,8e-4]],
#         [[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4],[8e-5,2e-4,4e-4,8e-4],[8e-5,2e-4,4e-4,8e-4]]]
# yticks=[[[np.log10(yticks[i][j][k]) for k in range(len(yticks[i][j]))] for j in range(nT)] for i in range(nmd)]
# yticklabels=[[[f'{10**yticks[i][j][k]:.0e}'.replace('e+0','e').replace('e-0','e-') for k in range(len(yticks[i][j]))] for j in range(nT)] for i in range(nmd)]
# yticks=np.log10(np.array([8e-5,2e-4,4e-4,8e-4]))
# yticklabels=[f'{10**yticks[i]:.0e}'.replace('e+0','e').replace('e-0','e-') for i in range(len(yticks))]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        for k in range(100):
            rgb=cmap(norm(k))
            ax.plot(np.arange(nch-1)+1,tau[i,j,k*100]*1e6,alpha=0.1)
        ax.plot(np.arange(nch-1)+1,tau_mean[i,j]*1e6,'k',linewidth=wth)
        #     ax.plot(np.arange(nch-1)+1,np.log10(tau[i,j,k*100]),alpha=0.1)
        # ax.plot(np.arange(nch-1)+1,np.log10(tau_mean[i,j]),'k',linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(45))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        # ax.set_yticks(yticks,yticklabels)
        if i==nmd-1:
            ax.set_xlabel('Mode\nat %sK'%int(Ts[j]),fontsize=size)
        if j==0:
            ax.set_ylabel('Relaxtion Time\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        ax.set_xscale('log')
        # ax.set_yscale('log')

ylim=[np.inf,-np.inf]
for ax in axs.flat:
    y1,y2=ax.get_ylim()
    ylim[0]=min(ylim[0],y1)
    ylim[1]=max(ylim[1],y2)
for ax in axs.flat:
    ax.set_ylim(ylim[0],ylim[1])

plt.tight_layout()
plt.savefig(f'{figname1}.png',format='png')
plt.savefig(f'{figname1}.pdf',format='pdf')
plt.show()

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
# yticks=np.log10(np.array([8e-5,2e-4,4e-4,8e-4]))
# yticklabels=[f'{10**yticks[i]:.0e}'.replace('e+0','e').replace('e-0','e-') for i in range(len(yticks))]

for i in range(nT):
    ax=axs[i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        # ax.plot(np.arange(nch-1)+1,np.log10(tau_mean[j,i]),color=rgb,linewidth=wth)
        ax.plot(np.arange(nch-1)+1,tau_mean[j,i]*1e6,color=rgb,linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(45))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_yticks(yticks,yticklabels)
    ax.set_xlabel('Mode\nat %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Relaxtion Time',fontsize=size)
    ax.set_xscale('log')
    # ax.set_yscale('log')

ylim=[np.inf,-np.inf]
for ax in axs.flat:
    y1,y2=ax.get_ylim()
    ylim[0]=min(ylim[0],y1)
    ylim[1]=max(ylim[1],y2)
for ax in axs.flat:
    ax.set_ylim(ylim[0],ylim[1])

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

plt.savefig(f'{figname1}_1.png',format='png')
plt.savefig(f'{figname1}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(15,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.88
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# yticks=np.log10(np.array([8e-5,2e-4,4e-4,8e-4]))
# yticklabels=[f'{10**yticks[i]:.0e}'.replace('e+0','e').replace('e-0','e-') for i in range(len(yticks))]

for i in range(nmd):
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        # ax.plot(np.arange(nch-1)+1,np.log10(tau_mean[i,j]),color=rgb,linewidth=wth)
        ax.plot(np.arange(nch-1)+1,tau_mean[i,j]*1e6,color=rgb,linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(45))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_yticks(yticks,yticklabels)
    ax.set_xlabel('Mode\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Relaxtion Time',fontsize=size)
    ax.set_xscale('log')
    # ax.set_yscale('log')

ylim=[np.inf,-np.inf]
for ax in axs.flat:
    y1,y2=ax.get_ylim()
    ylim[0]=min(ylim[0],y1)
    ylim[1]=max(ylim[1],y2)
for ax in axs.flat:
    ax.set_ylim(ylim[0],ylim[1])

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

plt.savefig(f'{figname1}_2.png',format='png')
plt.savefig(f'{figname1}_2.pdf',format='pdf')
plt.show()

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
        ax.plot(omg*1e-6,G1[j,i],color=rgb,linewidth=wth)
        ax.plot(omg*1e-6,G2[j,i],color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],'k',linewidth=wth,label='Elasticity')
    ax.plot([],[],'k--',linewidth=wth,label='Viscosity')
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
    ax.set_xlabel('Ang Freq\nat %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Modulus',fontsize=size)
        ax.legend(loc='best',fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

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

plt.savefig(f'{figname2}_1.png',format='png')
plt.savefig(f'{figname2}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(15,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.88
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

for i in range(nmd):
    ax=axs[i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(omg*1e-6,G1[i,j],color=rgb,linewidth=wth)
        ax.plot(omg*1e-6,G2[i,j],color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],'k',linewidth=wth,label='Elasticity')
    ax.plot([],[],'k--',linewidth=wth,label='Viscosity')
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
    ax.set_xlabel('Ang Freq\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Modulus\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        ax.legend(loc='best',fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

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

plt.savefig(f'{figname2}_2.png',format='png')
plt.savefig(f'{figname2}_2.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,2,figsize=(10,4))
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

ax=axs[0]  
spc=3
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.errorbar(np.arange(nT)*(nmd+1+spc)+i,X01_mean[i,:]*1e6,yerr=X01_std[i,:]*1e6,
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
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+1+spc)+(nmd-1)/2,np.array(Ts).astype(int))
ax.set_xlabel('Temperature (K)',fontsize=size)
ax.set_ylabel('0-Shear Compliance',fontsize=size)

ax=axs[1]  
spc=3
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.errorbar(np.arange(nT)*(nmd+1+spc)+i,X02_mean[i,:]*1e6,yerr=X02_std[i,:]*1e6,
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
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+1+spc)+(nmd-1)/2,np.array(Ts).astype(int))
ax.set_xlabel('Temperature (K)',fontsize=size)
ax.set_ylabel('0-Shear Viscosity',fontsize=size)

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

plt.savefig(f'{figname3}_1.png',format='png')
plt.savefig(f'{figname3}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,2,figsize=(10,4))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.83
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=axs[0]  
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+1+spc)+i,X01_mean[:,i]*1e6,yerr=X01_std[:,i]*1e6,
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
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+1+spc)+(nT-1)/2,np.around(10**rg_flu,2))
ax.set_xlabel('Structure Fluctuation',fontsize=size)
ax.set_ylabel('0-Shear Compliance',fontsize=size)

ax=axs[1]  
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+1+spc)+i,X02_mean[:,i]*1e6,yerr=X02_std[:,i]*1e6,
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
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+1+spc)+(nT-1)/2,np.around(10**rg_flu,2))
ax.set_xlabel('Structure Fluctuation',fontsize=size)
ax.set_ylabel('0-Shear Viscosity',fontsize=size)

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

plt.savefig(f'{figname3}_2.png',format='png')
plt.savefig(f'{figname3}_2.pdf',format='pdf')
plt.show()