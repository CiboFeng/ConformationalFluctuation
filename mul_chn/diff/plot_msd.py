"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
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
dir='msd'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
cut=[10,100]
dt=10e-6*10000 ### ns

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)
t=np.arange(nfr)*dt
d=np.zeros((nmd,nT,nfr,nch,3))
for i in range(nmd):
    for j in range(nT):
        d[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npy')

dxy=np.sum(d[:,:,:,:,:2],axis=-1)
dz=d[:,:,:,:,-1]
dm=np.mean(d,axis=-2)
dmxy=np.sum(dm[:,:,:,:2],axis=-1)
dmz=dm[:,:,:,-1]

dxy[:,:,-100:]=np.nan
dz[:,:,-100:]=np.nan
dmxy[:,:,-100:]=np.nan
dmz[:,:,-100:]=np.nan

"""Fit"""
def func(x,a,b):
    return a*x**b

Parxy=np.ones((nmd,nT,2))
Parz=np.ones((nmd,nT,2))
for i in range(nmd):
    for j in range(nT):
        par,cov=curve_fit(func,t[cut[0]:cut[1]],dmxy[i,j,cut[0]:cut[1]])
        Parxy[i,j]=par
        par,cov=curve_fit(func,t[cut[0]:cut[1]],dmz[i,j,cut[0]:cut[1]])
        Parz[i,j]=par

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(20,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nch)

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        for k in range(nch):
            rgb=cmap(norm(k))
            ax.plot(t,dxy[i,j,:,k],color=rgb,linewidth=wth,alpha=0.5)
            ax.plot(t,dz[i,j,:,k],color=rgb,linewidth=wth,alpha=0.5)
        ax.plot(t,dmxy[i,j,:],color='k',linewidth=wth)
        ax.plot(t,dmz[i,j,:],color=(0.5,0.5,0.5),linewidth=wth)
        ax.plot(t[cut[0]:cut[1]],3*Parxy[i,j,0]*t[cut[0]:cut[1]]**Parxy[i,j,1],color='k',linestyle=':',linewidth=wth,label='$\\langle{}\\mathrm{\\Delta}x^2+\\mathrm{\\Delta}y^2\\rangle{}=%.0f\\mathrm{\\Delta}t^{%.2f}$'%(Parxy[i,j,0],Parxy[i,j,1]))
        ax.plot(t[cut[0]:cut[1]],3*Parz[i,j,0]*t[cut[0]:cut[1]]**Parz[i,j,1],color=(0.5,0.5,0.5),linestyle=':',linewidth=wth,label='$\\langle{}\\mathrm{\\Delta}z^2\\rangle{}=%.0f\\mathrm{\\Delta}t^{%.2f}$'%(Parz[i,j,0],Parz[i,j,1]))
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        if i==nmd-1:
            ax.set_xlabel('Lag Time (ns) at %s K'%int(Ts[j]),fontsize=size)
        if j==0:
            ax.set_ylabel('MSD ($\\mathrm{\\AA}^2$) of Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_ylim(40,1500)
        ax.legend(loc='upper left',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(2,nT,figsize=(20,10))
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
    ax=axs[0,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(t,dmxy[j,i],color=rgb,linewidth=wth)
        ax.plot(t,dmz[j,i],color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='XY')
    ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='Z')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    if i==0:
        ax.set_ylabel('MSD ($\\mathrm{\\AA}^2$)',fontsize=size)
        ax.legend(loc='upper left',fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax=axs[1,i]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(t,dmxy[j,i]/t/6,color=rgb,linewidth=wth)
        ax.plot(t,dmz[j,i]/t/6,color=rgb,linestyle='--',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('Lag Time (ns) at %i K'%Ts[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Diffusion Coefficient ($\\mathrm{\\AA}^2$/ns)',fontsize=size)
    ax.set_xscale('log')
    # ax.set_yscale('log')

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

plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()
