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
dir='../../cont_life/cont_life'
figname='cont_life'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=35

"""Read Data and Calculate"""
T0=np.zeros((nmd,nT,nhi))
p_T0=np.zeros((nmd,nT,nhi))
T0_mean=np.zeros((nmd,nT))
T0_std=np.zeros((nmd,nT))
T0_stdl=np.zeros((nmd,nT))
T0_stdr=np.zeros((nmd,nT))
T1=np.zeros((nmd,nT,nhi))
p_T1=np.zeros((nmd,nT,nhi))
T1_mean=np.zeros((nmd,nT))
T1_std=np.zeros((nmd,nT))
T1_stdl=np.zeros((nmd,nT))
T1_stdr=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        T0[i,j]=q['T0']
        p_T0[i,j]=q['p_T0']
        T0_mean[i,j]=q['T0_mean']
        T0_std[i,j]=q['T0_std']
        T0_stdl[i,j]=q['T0_stdl']
        T0_stdr[i,j]=q['T0_stdr']
        T1[i,j]=q['T1']
        p_T1[i,j]=q['p_T1']
        T1_mean[i,j]=q['T1_mean']
        T1_std[i,j]=q['T1_std']
        T1_stdl[i,j]=q['T1_stdl']
        T1_stdr[i,j]=q['T1_stdr']

"""Plot"""
fig,axs=plt.subplots(2,nT+1,figsize=(25,10))
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
xtick=0.1
ytick=20
ctick=1
ann=[['(A)','(B)','(C)','(D)','(E)'],['(F)','(G)','(H)','(I)','(J)']]
clrs=['r','g','b']

for i in range(nT):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(T0[j,i,1:],p_T0[j,i,1:],color=clrs[j],linewidth=wth,label=Mdls[j]) ### Delete the sampling at lower limit.
        ax.plot(T0[j,i],[np.median(p_T0[j,i])]*nhi,linewidth=0.0)
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
    ax.set_xlabel('$\\tau_0$ (ns)',fontdict=font)
    ax.set_ylabel('$f(\\tau_0)~(\\mathrm{ns}^{-1})$',fontdict=font)
    if i==0:
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(T1[j,i,1:],p_T1[j,i,1:],color=clrs[j],linewidth=wth)
        ax.plot(T1[j,i],[np.median(p_T1[j,i])]*nhi,linewidth=0.0)
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
    ax.set_xlabel('$\\tau_1$ (ns) at %i K'%Ts[i],fontdict=font)
    ax.set_ylabel('$f(\\tau_1)~(\\mathrm{ns}^{-1})$',fontdict=font)

ylims=[axs[i,j].get_ylim() for j in range(nT) for i in range(2)]
for i in range(2):
    for j in range(nT):
        axs[i,j].set_ylim(np.min(ylims),np.max(ylims))

ax=axs[0,nT]  
ax.annotate(ann[0][nT],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,T0_mean[i,:],yerr=[(1-1/T0_stdl[i,:])*T0_mean[i,:],(T0_stdr[i,:]-1)*T0_mean[i,:]], ### This "standard deviation" is related to */, not +-.
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+spc)+(nmd-1)/2,np.array(Ts).astype(int),fontdict=font)
ax.set_ylabel('$\\tau_0$ (ns)',fontdict=font)

ax=axs[1,nT]  
ax.annotate(ann[1][nT],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,T1_mean[i,:],yerr=[(1-1/T1_stdl[i,:])*T1_mean[i,:],(T1_stdr[i,:]-1)*T1_mean[i,:]],
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+spc)+(nmd-1)/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$\\tau_1$ (ns)',fontdict=font)

ylims=[axs[i,-1].get_ylim() for i in range(2)]
for i in range(2):
    axs[i,-1].set_ylim(np.min(ylims),np.max(ylims))

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

plt.savefig(f'{figname}_1.png',format='png',dpi=50)
plt.savefig(f'{figname}_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nmd+1,figsize=(25,10))
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
xtick=0.1
ytick=20
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)']]

for i in range(nmd):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(T0[i,j,1:],p_T0[i,j,1:],color=cmap(j/(nT-1)),linewidth=wth,label='%i K'%Ts[j])
        ax.plot(T0[i,j],[np.median(p_T0[i,j])]*nhi,linewidth=0.0)
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
    ax.set_xlabel('$\\tau_0$ (ns)',fontdict=font)
    ax.set_ylabel('$f(\\tau_0)~(\\mathrm{ns}^{-1})$',fontdict=font)
    if i==0:
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(T1[i,j,1:],p_T1[i,j,1:],color=cmap(j/(nT-1)),linewidth=wth)
        ax.plot(T1[i,j],[np.median(p_T1[i,j])]*nhi,linewidth=0.0)
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
    ax.set_xlabel('$\\tau_1$ (ns) of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$f(\\tau_1)~(\\mathrm{ns}^{-1})$',fontdict=font)

ax=axs[0,nmd]  
ax.annotate(ann[0][nmd],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+spc)+i,T0_mean[:,i],yerr=[(1-1/T0_stdl[:,i])*T0_mean[:,i],(T0_stdr[:,i]-1)*T0_mean[:,i]],
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+spc)+(nT-1)/2,Mdls,fontdict=font)
ax.set_ylabel('$\\tau_0$ (ns)',fontdict=font)

ax=axs[1,nmd]  
ax.annotate(ann[1][nmd],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+spc)+i,T1_mean[:,i],yerr=[(1-1/T1_stdl[:,i])*T1_mean[:,i],(T1_stdr[:,i]-1)*T1_mean[:,i]],
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+spc)+(nT-1)/2,Mdls,fontdict=font)
ax.set_ylabel('$\\tau_1$ (ns)',fontdict=font)

ylims=[axs[i,-1].get_ylim() for i in range(2)]
for i in range(2):
    axs[i,-1].set_ylim(np.min(ylims),np.max(ylims))
    
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

plt.savefig(f'{figname}_2.png',format='png',dpi=50)
plt.savefig(f'{figname}_2.pdf',format='pdf')
# plt.show()