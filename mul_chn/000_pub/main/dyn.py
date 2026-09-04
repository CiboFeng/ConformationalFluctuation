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
dir1='../../diff/msd_in'
file2='../../autocor/proc_autocor_dist_in.npz'
dir3='../../cont_life/cont_life'
dir4='../../strs/grf_mod'
figname='dyn'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
dt=10e-6*10000 ### ns
nat=163
nch=100
nhi=35
nomg=100

"""Diffusion"""
t=np.arange(nfr)*dt
d=np.zeros((nmd,nT,nfr,nch,3))
for i in range(nmd):
    for j in range(nT):
        d[i,j]=np.load(f'{dir1}/{mdls[i]}_{Ts[j]}.npy')
d[d==0.0]=np.nan
d[:,:,0]=0.0

dxy=np.sum(d[:,:,:,:,:2],axis=-1)
dz=d[:,:,:,:,-1]
dm=np.nanmean(d,axis=-2)
dmxy=np.sum(dm[:,:,:,:2],axis=-1)
dmz=dm[:,:,:,-1]

def lastnan(X,last):
    num_idx=np.where(~np.isnan(X))[0]
    last_idx=num_idx[-last:] if len(num_idx)>=last else num_idx
    X[last_idx]=np.nan
for i in range(nmd):
    for j in range(nT):
        for k in range(nch):
            lastnan(dxy[i,j,:,k],100)
            lastnan(dz[i,j,:,k],100)
        lastnan(dmxy[i,j],100)
        lastnan(dmz[i,j],100)

"""Distance Relaxation Time"""
t=np.arange(nfr)*dt
q=np.load(file2)
ac_mean=q['ac_mean']
t_relx_mean=q['t_relx_mean']
t_relxe_mean=q['t_relxe_mean']

"""Contact Life Time"""
T0=np.zeros((nmd,nT,nhi))
p_T0=np.zeros((nmd,nT,nhi))
T0_mean=np.zeros((nmd,nT))
T0_std=np.zeros((nmd,nT))
T1=np.zeros((nmd,nT,nhi))
p_T1=np.zeros((nmd,nT,nhi))
T1_mean=np.zeros((nmd,nT))
T1_std=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir3}/{mdls[i]}_{Ts[j]}.npz')
        T0[i,j]=q['T0']
        p_T0[i,j]=q['p_T0']
        T0_mean[i,j]=q['T0_mean']
        T0_std[i,j]=q['T0_std']
        T1[i,j]=q['T1']
        p_T1[i,j]=q['p_T1']
        T1_mean[i,j]=q['T1_mean']
        T1_std[i,j]=q['T1_std']

"""Stress"""
tau=np.zeros((nmd,nT,nfr,nch-1))
for i in range(nmd):
    for j in range(nT):
        tau[i,j]=np.load(f'{dir4}/{mdls[i]}_{Ts[j]}.npy')
tau_mean=np.mean(tau,axis=2)
# omg=np.logspace(3,7,nomg)
omg=np.linspace(0,2e6,nomg)
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
fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3)
Gs=[gs[0].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3),
    gs[1].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3),
    gs[2].subgridspec(nrows=3,ncols=1,height_ratios=[1,1,1],hspace=0.4)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].nrows)] for i in range(len(Gs))]
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=-0.1
cbarwth=0.0
cbarspcwth=0.0
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
lenbar=8
xtick=100
ytick=100
ctick=1
clrs=['r','g','b']
ann=[['(A)','(B)'],['(C)','(D)'],['(E)','(F)','(G)']]
clrs_=['#843204','#EF00C2']

ax=axs[0][0]
ax.annotate(ann[0][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t,dmz[i,0]/t/6,color=clrs[i],linewidth=wth,label=Mdls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$\\Delta{}t~(\\mathrm{ns})$',fontdict=font)
ax.set_ylabel('$D~(\\mathrm{\\AA^2/ns})$',fontdict=font)
ax.set_xscale('log')
# ax.set_yscale('log')
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ax=axs[0][1]
ax.annotate(ann[0][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(np.arange(nat),t_relxe_mean[i,0],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$|i-j|$',fontdict=font)
ax.set_ylabel('$\\tau(d_{ij})$ (ns)',fontdict=font)

ax=axs[1][0]
ax.annotate(ann[1][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(T1[i,0,1:],p_T1[i,0,1:],color=clrs[i],linewidth=wth)
    ax.plot(T1[i,0],[np.median(p_T1[i,0])]*nhi,linewidth=0.0)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_minor_locator(AutoMinorLocator(0))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([1e0,1e2],['$10^{0}$','$10^{2}$'],fontdict=font)
ax.set_xlabel('$\\tau_1$ (ns)',fontdict=font)
ax.set_ylabel('$f(\\tau_1)~(\\mathrm{ns^{-1}})$',fontdict=font)

ax=axs[1][1]
ax.annotate(ann[1][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    # ax.errorbar(i,T1_mean[i,0],yerr=T1_std[i,0],
    #                  linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
    #                  ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
    ax.bar(i,T1_mean[i,0],width=0.3,color=clrs[i],edgecolor='k',linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nmd),Mdls,fontdict=font)
ax.set_ylabel('$\\langle\\tau_1\\rangle$ (ns)',fontdict=font)

ax=axs[2][0]
ax.annotate(ann[2][0],xy=(0.0,1.075),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(np.arange(nch-1)+1,tau_mean[i,0]*1e6,color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_yticks(yticks,yticklabels,fontdict=font)
ax.set_xlabel('$n_\\mathrm{mod}$',fontdict=font)
ax.set_ylabel('$\\langle\\tau_n\\rangle$',fontdict=font)
ax.set_xscale('log')
# ax.set_yscale('log')

ax=axs[2][1]
ax.annotate(ann[2][1],xy=(0.0,1.075),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(omg*1e-6,G1[i,0],color=clrs[i],linewidth=wth)
    ax.plot(omg*1e-6,G2[i,0],color=clrs[i],linestyle='--',linewidth=wth)
ax.plot([],[],'k',linewidth=wth,label='$G^\\prime$')
ax.plot([],[],'k--',linewidth=wth,label='$G^{\\prime\\prime}$')
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.4))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$\\omega$',fontdict=font)
ax.set_ylabel('$G^\\prime~&~G^{\\prime\\prime}$',fontdict=font)
ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
# ax.set_xscale('log')
# ax.set_yscale('log')

ax=axs[2][2]
ax.annotate(ann[2][2],xy=(0.0,1.075),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
ax2=ax.twinx()
spc=0.3
y1min=1.1*np.min(X01_mean[:,0]-X01_std[:,0])-0.1*np.max(X01_mean[:,0]+X01_std[:,0])
y2min=1.1*np.min(X02_mean[:,0]-X02_std[:,0])-0.1*np.max(X02_mean[:,0]+X02_std[:,0])
X01_mean[:,0]-=y1min
X02_mean[:,0]-=y2min
for i in range(nmd):
    ax.bar(i-spc/2,X01_mean[i,0]*1e6,width=spc,color=clrs[i],edgecolor=clrs_[0],linewidth=2*wth,alpha=0.5)
    ax.bar(i-spc/2,X01_mean[i,0]*1e6,width=spc,color='none',edgecolor=clrs_[0],linewidth=2*wth)
    ax.errorbar(i-spc/2,X01_mean[i,0]*1e6,yerr=X01_std[i,0]*1e6,
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
    ax2.bar(i+spc/2,X02_mean[i,0]*1e6,width=spc,color=clrs[i],edgecolor=clrs_[1],linewidth=2*wth,alpha=0.5)
    ax2.bar(i+spc/2,X02_mean[i,0]*1e6,width=spc,color='none',edgecolor=clrs_[1],linewidth=2*wth)
    ax2.errorbar(i+spc/2,X02_mean[i,0]*1e6,yerr=X02_std[i,0]*1e6,
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'],color=clrs_[0],labelcolor=clrs_[0])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'],color=clrs_[0],labelcolor=clrs_[0])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.spines['left'].set_color(clrs_[0])
ax.yaxis.label.set_color(clrs_[0])
ax.set_xticks(np.arange(nmd),Mdls,fontdict=font)
yticks=np.array([3400,3420,3440])
ax.set_yticks(yticks-y1min*1e6,yticks,fontdict=font)
ax.set_ylabel('$g_\\mathrm{e0}$',fontdict=font)
ax2.minorticks_on()
ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'],color=clrs_[1],labelcolor=clrs_[1])
ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'],color=clrs_[1],labelcolor=clrs_[1])
ax2.xaxis.set_minor_locator(AutoMinorLocator(0))
ax2.yaxis.set_major_locator(MultipleLocator(20))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(wth)
ax2.spines['right'].set_color(clrs_[1])
ax2.yaxis.label.set_color(clrs_[1])
yticks=np.array([300,320,340,360])
ax2.set_yticks(yticks-y2min*1e6,yticks,fontdict=font)
ax2.set_ylabel('$g_\\mathrm{v0}$',fontdict=font)

fig.canvas.draw()
tight_bbox=fig.get_tightbbox(fig.canvas.get_renderer())
alllft=tight_bbox.x0
allrig=tight_bbox.x1
allbot=tight_bbox.y0
alltop=tight_bbox.y1
figwth,fighei=fig.get_size_inches()
boxlft=np.inf
boxrig=-np.inf
boxbot=np.inf
boxtop=-np.inf
for ax in fig.axes:    
    pos=ax.get_position()
    left=pos.x0*figwth
    right=left+pos.width*figwth
    bottom=pos.y0*fighei
    top=bottom+pos.height*fighei
    boxlft=min(boxlft,left)
    boxrig=max(boxrig,right)
    boxbot=min(boxbot,bottom)
    boxtop=max(boxtop,top)
outlft=boxlft-alllft
outrig=allrig-boxrig
outbot=boxbot-allbot
outtop=alltop-boxtop

plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')