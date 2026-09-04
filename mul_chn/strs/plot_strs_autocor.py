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
dir='strs_autocor'
figname1='strs_autocor'
figname2='shear_mod'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=1000
dt=10 ### fs

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

ac=np.zeros((nmd,nT,6,nfr))
freq=np.zeros((nmd,nT,round(nfr/2)))
A=np.zeros((nmd,nT,6,round(nfr/2)))
phi=np.zeros((nmd,nT,6,round(nfr/2)))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        ac[i,j]=q['ac']
        freq[i,j]=q['freq']
        A[i,j]=q['A']
        phi[i,j]=q['phi']
G1=A*np.cos(phi)
G2=A*np.sin(phi)
G1_ave=np.mean(G1[:,:,3:],axis=2)
G2_ave=np.mean(G2[:,:,3:],axis=2)

"""Plot"""
fig,axs=plt.subplots(6,nT,figsize=(15,15))
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
ylabels=['$P_{xx}$','$P_{yy}$','$P_{zz}$','$P_{xy}$','$P_{yz}$','$P_{zx}$']

for i in range(6):
    for j in range(nT):
        ax=axs[i,j]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(np.arange(1,len(ac[k,j,i]))*dt,ac[k,j,i,1:],color=rgb,linewidth=wth)
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
        if i==5:
            ax.set_xlabel('Lag Time (fs)\nat %i K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Autocorrelation\nof %s'%ylabels[i],fontsize=size)
        ax.set_xscale('log')

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

plt.savefig(f'{figname1}_1.png',format='png')
plt.savefig(f'{figname1}_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(6,nmd,figsize=(15,15))
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
ylabels=['$P_{xx}$','$P_{yy}$','$P_{zz}$','$P_{xy}$','$P_{yz}$','$P_{zx}$']

for i in range(6):
    for j in range(nmd):
        ax=axs[i,j]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(np.arange(1,len(ac[j,k,i]))*dt,ac[j,k,i,1:],color=rgb,linewidth=wth)
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
        if i==5:
            ax.set_xlabel('Lag Time (fs)\nof Flu.=%.2f'%10**rg_flu[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Autocorrelation\nof %s'%ylabels[i],fontsize=size)
        ax.set_xscale('log')

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

plt.savefig(f'{figname1}_2.png',format='png')
plt.savefig(f'{figname1}_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(6,nT,figsize=(15,15))
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
ylabels=['xy','yz','zx']

for i in range(3):
    for j in range(nT):
        ax=axs[2*i,j]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(freq[k,j],G1[k,j,i+3],color=rgb,linewidth=wth)
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
        ymin=np.min(G1[:,:,3:])
        ymax=np.max(G1[:,:,3:])
        ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
        if j==0:
            ax.set_ylabel('Elastic Modulus\nof %s'%ylabels[i],fontsize=size)
        ax.set_xscale('log')

for i in range(3):
    for j in range(nT):
        ax=axs[2*i+1,j]
        for k in range(nmd):
            rgb=cmap(k/(nmd-1))
            ax.plot(freq[k,j],G2[k,j,i+3],color=rgb,linewidth=wth)
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
        ymin=np.min(G2[:,:,3:])
        ymax=np.max(G2[:,:,3:])
        ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
        if i==2:
            ax.set_xlabel('Frequency (fs$^{-1}$)\nat %i K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Viscous Modulus\nof %s'%ylabels[i-3],fontsize=size)
        ax.set_xscale('log')

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

plt.savefig(f'{figname2}_1.png',format='png')
plt.savefig(f'{figname2}_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(6,nmd,figsize=(15,15))
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
ylabels=['xy','yz','zx']

for i in range(3):
    for j in range(nmd):
        ax=axs[2*i,j]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(freq[j,k],G1[j,k,i+3],color=rgb,linewidth=wth)
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
        ymin=np.min(G1[:,:,3:])
        ymax=np.max(G1[:,:,3:])
        ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
        if j==0:
            ax.set_ylabel('Elastic Modulus\nof %s'%ylabels[i],fontsize=size)
        ax.set_xscale('log')

for i in range(3):
    for j in range(nmd):
        ax=axs[2*i+1,j]
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(freq[j,k],G2[j,k,i+3],color=rgb,linewidth=wth)
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
        ymin=np.min(G2[:,:,3:])
        ymax=np.max(G2[:,:,3:])
        ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
        if i==2:
            ax.set_xlabel('Frequency (fs$^{-1}$)\nof Flu.=%.2f'%10**rg_flu[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Viscous Modulus\nof %s'%ylabels[i-3],fontsize=size)
        ax.set_xscale('log')

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

plt.savefig(f'{figname2}_2.png',format='png')
plt.savefig(f'{figname2}_2.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nT,figsize=(15,6))
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
ylabels=['xy','yz','zx']

for i in range(nT):
    ax=axs[0,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(freq[j,i],G1_ave[j,i],color=rgb,linewidth=wth)
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
    ymin=np.min(G1_ave)
    ymax=np.max(G1_ave)
    ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
    if i==0:
        ax.set_ylabel('Elastic Modulus',fontsize=size)
    ax.set_xscale('log')

for i in range(nT):
    ax=axs[1,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(freq[j,i],G2_ave[j,i],color=rgb,linewidth=wth)
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
    ymin=np.min(G2_ave)
    ymax=np.max(G2_ave)
    ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
    ax.set_xlabel('Frequency (fs$^{-1}$)\nat %i K'%Ts[j],fontsize=size)
    if i==0:
        ax.set_ylabel('Viscous Modulus',fontsize=size)
    ax.set_xscale('log')

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

plt.savefig(f'{figname2}_ave_1.png',format='png')
plt.savefig(f'{figname2}_ave_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(2,nmd,figsize=(15,6))
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
ylabels=['xy','yz','zx']

for i in range(nmd):
    ax=axs[0,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(freq[i,j],G1_ave[i,j],color=rgb,linewidth=wth)
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
    ymin=np.min(G1_ave)
    ymax=np.max(G1_ave)
    ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
    if i==0:
        ax.set_ylabel('Elastic Modulus',fontsize=size)
    ax.set_xscale('log')

for i in range(nmd):
    ax=axs[1,i]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(freq[i,j],G2_ave[i,j],color=rgb,linewidth=wth)
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
    ymin=np.min(G2_ave)
    ymax=np.max(G2_ave)
    ax.set_ylim(1.05*ymin-0.05*ymax,-0.05*ymin+1.05*ymax)
    ax.set_xlabel('Frequency (fs$^{-1}$)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Viscous Modulus',fontsize=size)
    ax.set_xscale('log')

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

plt.savefig(f'{figname2}_ave_2.png',format='png')
plt.savefig(f'{figname2}_ave_2.pdf',format='pdf')
plt.show()