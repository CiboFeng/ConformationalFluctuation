"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
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
dir='autocor_cont_seq'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nfr=10000
nat=163

"""Read Data"""
t=np.zeros((nmd,nfr,nat))
ac=np.zeros((nmd,nfr,nat))
for i in range(nmd):
    q=pyw(f'{dir}/{mdls[i]}_dyn.pyw','Autocorrelation')
    for j in range(nfr):
        for k in range(nat):
            t[i,j,k]=q[0,j*nat+k]
            ac[i,j,k]=q[2,j*nat+k]

"""Interpolation"""
ac_tgt=np.arange(0.1,1.0,0.1)
t_relx=np.zeros((nmd,nat,len(ac_tgt)))
t_relxe=np.zeros((nmd,nat))
for i in range(nmd):
    for j in range(nat):
        for k in range(len(ac_tgt)):
            # t_relx[i,j,k]=np.interp(ac_tgt[k],ac[i,:,j],t[i,:,j])
            for l in range(nfr-1):
                if ac[i,l,j]>ac_tgt[k] and ac[i,l+1,j]<=ac_tgt[k]:
                    t_relx[i,j,k]=((ac[i,l+1,j]-ac_tgt[k])*t[i,l,j]+(ac_tgt[k]-ac[i,l,j])*t[i,l+1,j])/(ac[i,l+1,j]-ac[i,l,j])
                    break
        for l in range(nfr-1):
            if ac[i,l,j]>1/np.e and ac[i,l+1,j]<=1/np.e:
                t_relxe[i,j]=((ac[i,l+1,j]-1/np.e)*t[i,l,j]+(1/np.e-ac[i,l,j])*t[i,l+1,j])/(ac[i,l+1,j]-ac[i,l,j])
                break
t_relxe[:,0]=np.nan

"""Fit"""
# def func(x,a,b):
#     return a*x**b

# Par=np.ones((nc,2))
# for i in range(nc):
#     par,cov=curve_fit(func,t[i,0,:cut],dmm[i,:cut])
#     Par[i]=par

"""Plot"""
fig,axs=plt.subplots(1,3,figsize=(15,4))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)

for i in range(nmd):
    ax=axs[i]
    for j in range(nat):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(t[i,:,j],ac[i,:,j],color=rgb,linewidth=wth)
    ax.plot([],[],label=labels[i])
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
    ax.set_xlabel('Lag Time (ns)',fontsize=size)
    if i==0:
        ax.set_ylabel(f'Autocorrelation',fontsize=size)
    ax.set_xscale('log')
    ax.set_ylim(-0.05,1.05)
    ax.legend(loc='upper right',fontsize=size,handlelength=0.0,handletextpad=0.0)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_1.png',format='png')
plt.savefig(f'{figname}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,3,figsize=(15,4))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.LogNorm(vmin=np.min(t[:,1:101,:]),vmax=np.max(t[:,1:101,:]))

for i in range(nmd):
    ax=axs[i]
    for j in range(len(t[i,1:101,0])):
        normc=norm(t[i,j,0])
        rgb=cmap(normc)
        ax.plot(ac[i,j,:],color=rgb,linewidth=wth)
    ax.plot([],[],label=labels[i])
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
    ax.set_xlabel('Sequence Distance',fontsize=size)
    if i==0:
        ax.set_ylabel(f'Autocorrelation',fontsize=size)
    # ax.set_xscale('log')
    ax.set_ylim(-0.05,1.05)
    ax.legend(loc='upper right',fontsize=size,handlelength=0.0,handletextpad=0.0)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Lag Time (ns)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(1,3,figsize=(15,4))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=np.min(ac_tgt),vmax=np.max(ac_tgt))

for i in range(nmd):
    ax=axs[i]
    for j in range(len(ac_tgt)):
        normc=norm(ac_tgt[j])
        rgb=cmap(normc)
        ax.plot(t_relx[i,:,j],color=rgb,linewidth=wth)
    ax.plot([],[],label=labels[i])
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
    ax.set_xlabel('Sequence Distance',fontsize=size)
    if i==0:
        ax.set_ylabel('Relaxation Time (ns)',fontsize=size)
    ax.set_yscale('log')
    ax.set_ylim(8e-3,2.5e0)
    ax.legend(loc='lower center',fontsize=size,handlelength=0.0,handletextpad=0.0)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Autocorrelation',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_3.png',format='png')
plt.savefig(f'{figname}_3.pdf',format='pdf')
plt.show()

fig=plt.figure(figsize=(6,4))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
colors=['r','g','b']
dseq=np.arange(nat)

ax=plt.subplot(111)
for i in range(nmd):
    ax.plot(dseq,t_relxe[i],color=colors[i],linewidth=wth,label=labels[i])
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
ax.set_xlabel('Sequence Distance',fontsize=size)
ax.set_ylabel('Relaxation Time (ns)',fontsize=size)
# ax.set_yscale('log')
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_4.png',format='png')
plt.savefig(f'{figname}_4.pdf',format='pdf')
plt.show()