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
dir='../../strs/grf_mod'
figname='grf_mod'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
nomg=100

"""Read Data"""
tau=np.zeros((nmd,nT,nfr,nch-1))
for i in range(nmd):
    for j in range(nT):
        tau[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npy')
tau_mean=np.mean(tau,axis=2)
tau_std=np.std(tau,axis=2)
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
fig,axs=plt.subplots(nmd,nT,figsize=(25,15))
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
xtick=0
ytick=1
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=100)
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        # for k in range(100):
        #     ax.plot(np.arange(nch-1)+1,tau[i,j,k*100]*1e6,color=cmap(norm(k)),alpha=0.1)
        ax.fill_between(np.arange(nch-1)+1,tau_mean[i,j]*1e6-tau_std[i,j]*1e6,tau_mean[i,j]*1e6+tau_std[i,j]*1e6,color='k',linewidth=0,alpha=0.2)
        ax.plot(np.arange(nch-1)+1,tau_mean[i,j]*1e6,'k',linewidth=wth)
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_minor_locator(AutoMinorLocator(0))
        ax.yaxis.set_major_locator(MultipleLocator(ytick))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        if i==nmd-1:
            ax.set_xlabel('$n$ at %i K'%Ts[j],fontdict=font)
        if j==0:
            ax.set_ylabel('$\\tau_n$ of %s'%Mdls[i],fontdict=font)
        ax.set_xscale('log')

ylim=[np.inf,-np.inf]
for ax in axs.flat:
    y1,y2=ax.get_ylim()
    ylim[0]=min(ylim[0],y1)
    ylim[1]=max(ylim[1],y2)
for ax in axs.flat:
    ax.set_ylim(ylim[0],ylim[1])

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

axs_2d=np.atleast_2d(axs)
nrows,ncols=np.shape(axs_2d)
outlfts=np.zeros((nrows,ncols))
outrigs=np.zeros((nrows,ncols))
outbots=np.zeros((nrows,ncols))
outtops=np.zeros((nrows,ncols))
dpi=fig.dpi 
for i in range(nrows):
    for j in range(ncols):
        ax=axs_2d[i,j]
        ref_bounds=ax.get_position().bounds
        alllft,allrig=np.inf,-np.inf
        allbot,alltop=np.inf,-np.inf
        for a in fig.axes:
            if a.get_position().bounds==ref_bounds:
                tight_bbox=a.get_tightbbox(fig.canvas.get_renderer())   # 点坐标
                alllft=min(alllft,tight_bbox.x0/dpi)
                allrig=max(allrig,tight_bbox.x1/dpi)
                allbot=min(allbot,tight_bbox.y0/dpi)
                alltop=max(alltop,tight_bbox.y1/dpi)
        boxlft=ref_bounds[0]*figwth
        boxrig=(ref_bounds[0]+ref_bounds[2])*figwth
        boxbot=ref_bounds[1]*fighei
        boxtop=(ref_bounds[1]+ref_bounds[3])*fighei
        outlfts[i,j]=boxlft-alllft
        outrigs[i,j]=allrig-boxrig
        outbots[i,j]=boxbot-allbot
        outtops[i,j]=alltop-boxtop

if ncols>1:
    spcwth+=np.max(outrigs[:,:-1]+outlfts[:,1:])
if nrows>1:
    spchei+=np.max(outbots[:-1,:]+outtops[1:,:])
plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei,
                    hspace=nrows*spchei/(fighei-margbot-outbot-margtop-outtop-(nrows-1)*spchei),
                    wspace=ncols*spcwth/(figwth-marglft-outlft-margrig-outrig-cbarwth-cbarspcwth-(ncols-1)*spcwth))

plt.savefig(f'{figname}_relx_time.png',format='png',dpi=50)
plt.savefig(f'{figname}_relx_time.pdf',format='pdf')
# plt.show()

fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=3,ncols=1,height_ratios=[1,1,1],hspace=0.4)
Gs=[gs[0].subgridspec(nrows=1,ncols=nT,width_ratios=[1]*nT,wspace=0.3),
    gs[1].subgridspec(nrows=1,ncols=nT,width_ratios=[1]*nT,wspace=0.3),
    gs[2].subgridspec(nrows=1,ncols=2,width_ratios=[1,1],wspace=0.132)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
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
xtick=[[],[1]*nT,[]]
ytick=[[1]*nT,[0.5]*nT,[500,100]]
ctick=1
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)']]
clrs=['r','g','b']

for i in range(nT):
    ax=axs[0][i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(np.arange(nch-1)+1,tau_mean[j,i]*1e6,color=clrs[j],linewidth=wth,label=Mdls[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[0][i]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_yticks(yticks,yticklabels,fontdict=font)
    ax.set_xlabel('$n$ at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$\\langle\\tau_n\\rangle$',fontdict=font)
        leg=ax.legend(loc='upper right',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    # ax.set_yscale('log')

    ax=axs[1][i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(omg*1e-6,G1[j,i],color=clrs[j],linewidth=wth)
        ax.plot(omg*1e-6,G2[j,i],color=clrs[j],linestyle='--',linewidth=wth)
    ax.plot([],[],'k',linewidth=wth,label='$G^\\prime$')
    ax.plot([],[],'k--',linewidth=wth,label='$G^{\\prime\\prime}$')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(xtick[1][i]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\omega$ at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$G^\\prime(\\omega)~&~G^{\\prime\\prime}(\\omega)$',fontdict=font)
        ax.legend(loc='lower right',prop={'family':font['family'],'size':font['size']})

for i in range(2):
    ylims=[axs[i][j].get_ylim() for j in range(nT)]
    for j in range(nT):
        axs[i][j].set_ylim(np.min(ylims),np.max(ylims))

ax=axs[2][0] 
ax.annotate(ann[2][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right') 
spc=3
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,X01_mean[i,:]*1e6,yerr=X01_std[i,:]*1e6,
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(ytick[2][0]))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+spc)+(nmd-1)/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$g_\\mathrm{e0}$',fontdict=font)

ax=axs[2][1]  
ax.annotate(ann[2][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,X02_mean[i,:]*1e6,yerr=X02_std[i,:]*1e6,
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(ytick[2][1]))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+spc)+(nmd-1)/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$g_\\mathrm{v0}$',fontdict=font)

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

plt.savefig(f'{figname}_1.png',format='png',dpi=50)
plt.savefig(f'{figname}_1.pdf',format='pdf')
# plt.show()

fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=3,ncols=1,height_ratios=[1,1,1],hspace=0.4)
Gs=[gs[0].subgridspec(nrows=1,ncols=nmd,width_ratios=[1]*nmd,wspace=0.15),
    gs[1].subgridspec(nrows=1,ncols=nmd,width_ratios=[1]*nmd,wspace=0.15),
    gs[2].subgridspec(nrows=1,ncols=2,width_ratios=[1,1],wspace=0.15)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
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
xtick=[[],[1]*nmd,[]]
ytick=[[1]*nmd,[0.5]*nmd,[500,100]]
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)'],['(G)','(H)']]

for i in range(nmd):
    ax=axs[0][i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(np.arange(nch-1)+1,tau_mean[i,j]*1e6,color=cmap(j/(nT-1)),linewidth=wth,label='%i K'%Ts[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[0][i]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_yticks(yticks,yticklabels,fontdict=font)
    ax.set_xlabel('$n$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$\\langle\\tau_n\\rangle$',fontdict=font)
        leg=ax.legend(loc='upper right',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    # ax.set_yscale('log')

    ax=axs[1][i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(omg*1e-6,G1[i,j],color=rgb,linewidth=wth)
        ax.plot(omg*1e-6,G2[i,j],color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],'k',linewidth=wth,label='$G^\\prime$')
    ax.plot([],[],'k--',linewidth=wth,label='$G^{\\prime\\prime}$')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(xtick[1][i]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\omega$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$G^\\prime(\\omega)~&~G^{\\prime\\prime}(\\omega)$',fontdict=font)
        ax.legend(loc='lower right',prop={'family':font['family'],'size':font['size']})

for i in range(2):
    ylims=[axs[i][j].get_ylim() for j in range(nmd)]
    for j in range(nmd):
        axs[i][j].set_ylim(np.min(ylims),np.max(ylims))

ax=axs[2][0] 
ax.annotate(ann[2][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right') 
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+spc)+i,X01_mean[:,i]*1e6,yerr=X01_std[:,i]*1e6,
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(ytick[2][0]))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+spc)+(nT-1)/2,Mdls,fontdict=font)
ax.set_ylabel('$g_\\mathrm{e0}$',fontdict=font)

ax=axs[2][1]  
ax.annotate(ann[2][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+spc)+i,X02_mean[:,i]*1e6,yerr=X02_std[:,i]*1e6,
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(ytick[2][1]))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+spc)+(nT-1)/2,Mdls,fontdict=font)
ax.set_ylabel('$g_\\mathrm{v0}$',fontdict=font)

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

plt.savefig(f'{figname}_2.png',format='png',dpi=50)
plt.savefig(f'{figname}_2.pdf',format='pdf')
# plt.show()