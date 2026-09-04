"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
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
dir='autocor_dist'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nfr=10000
nat=163
dt=10e-6*10000 ### ns

"""Read Data"""
t=np.arange(nfr)*dt
ac=np.zeros((nmd,nfr,nat,nat))
for i in range(nmd):
    ac[i]=np.load(f'{dir}/{mdls[i]}_dyn.npy')

"""Interpolation"""
ac_tgt=np.arange(0.1,1.0,0.1)
t_relx=np.zeros((nmd,nat,nat,len(ac_tgt)))
t_relxe=np.zeros((nmd,nat,nat))
for i in range(nmd):
    for j in range(nat):
        for k in range(j):
            for l in range(len(ac_tgt)):
                # t_relx[i,j,l]=np.interp(ac_tgt[l],ac[i,:,j],t[i,:,j])
                for m in range(nfr-1):
                    if ac[i,m,j,k]>ac_tgt[l] and ac[i,m+1,j,k]<=ac_tgt[l]:
                        t_relx[i,j,k,l]=((ac[i,m+1,j,k]-ac_tgt[l])*t[m]+(ac_tgt[l]-ac[i,m,j,k])*t[m+1])/(ac[i,m+1,j,k]-ac[i,m,j,k])
                        t_relx[i,k,j,l]=t_relx[i,j,k,l]
                        break
            for m in range(nfr-1):
                if ac[i,m,j,k]>1/np.e and ac[i,m+1,j,k]<=1/np.e:
                    t_relxe[i,j,k]=((ac[i,m+1,j,k]-1/np.e)*t[m]+(1/np.e-ac[i,m,j,k])*t[m+1])/(ac[i,m+1,j,k]-ac[i,m,j,k])
                    t_relxe[i,k,j]=t_relxe[i,j,k]
                    break

ac_mean=np.zeros((nmd,nfr,nat))
t_relx_mean=np.zeros((nmd,nat,len(ac_tgt)))
t_relxe_mean=np.zeros((nmd,nat))
for i in range(nat):
    for j in range(nat-i):
        ac_mean[:,:,i]+=ac[:,:,j,i+j]
        t_relx_mean[:,i]+=t_relx[:,j,i+j]
        t_relxe_mean[:,i]+=t_relxe[:,j,i+j]
    ac_mean[:,:,i]/=nat-i
    t_relx_mean[:,i]/=nat-i
    t_relxe_mean[:,i]/=nat-i
t_relxe_mean[:,0]=np.nan

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
npair=round(nat/2)
dseq=np.arange(1,nat,round(nat/npair))
IJ=[]
for i in range(len(dseq)):
    I=random.randint(0,nat-1-dseq[i])
    J=I+dseq[i]
    IJ+=[(I,J)]

for i in range(nmd):
    ax=axs[i]
    # for j,k in IJ:
    #     normc=norm(abs(j-k))
    #     rgb=cmap(normc)
    #     ax.plot(t,ac[i,:,j,k],color=rgb,linewidth=wth)
    for j in range(nat):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(t,ac_mean[i,:,j],color=rgb,linewidth=wth)
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
norm=colors.LogNorm(vmin=np.min(t[1:101]),vmax=np.max(t[1:101]))
dseq=np.arange(nat)
IJ=[]
for i in range(nat):
    I=random.randint(0,nat-1-dseq[i])
    J=I+dseq[i]
    IJ+=[[I,J]]
IJ=tuple(np.array(IJ).T)

for i in range(nmd):
    ax=axs[i]
    for j in range(len(t[1:101])):
        normc=norm(t[j])
        rgb=cmap(normc)
        # ax.plot(ac[i,j][IJ],color=rgb,linewidth=wth)
        ax.plot(ac_mean[i,j],color=rgb,linewidth=wth)
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
dseq=np.arange(nat)
IJ=[]
for i in range(nat):
    I=random.randint(0,nat-1-dseq[i])
    J=I+dseq[i]
    IJ+=[[I,J]]
IJ=tuple(np.array(IJ).T)

for i in range(nmd):
    ax=axs[i]
    for j in range(len(ac_tgt)):
        normc=norm(ac_tgt[j])
        rgb=cmap(normc)
        # ax.plot(dseq,t_relx[i][IJ][:,j],color=rgb,linewidth=wth)
        ax.plot(dseq,t_relx_mean[i,:,j],color=rgb,linewidth=wth)
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
    ax.set_ylim(7e-3,6e0)
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
IJ=[]
for i in range(nat):
    I=random.randint(0,nat-1-dseq[i])
    J=I+dseq[i]
    IJ+=[[I,J]]
IJ=tuple(np.array(IJ).T)

ax=plt.subplot(111)
for i in range(nmd):
    # ax.plot(dseq,t_relxe[i][IJ],color=colors[i],linewidth=wth,label=labels[i])
    ax.plot(dseq,t_relxe_mean[i],color=colors[i],linewidth=wth,label=labels[i])
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