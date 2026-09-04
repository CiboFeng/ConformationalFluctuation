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
dir='../../b22/pmf_rep'
figname='b22'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
nhi=170
nrep=3
kB=0.001982923700
T=300.0

"""Read Data"""
d=np.zeros((nmd,nhi,nrep))
F=np.zeros((nmd,nhi,nrep))
for i in range(nmd):
    for j in range(nrep):
        f=open(f'{dir}/{mdls[i]}_{j}.dat','r')
        ls=f.readlines()
        k=0
        for l in range(len(ls)):
            l_s=ls[l].strip('\n').split()
            if len(l_s)==5 and l_s[0][0]!='#':
                d[i,k,j]=float(l_s[0])
                F[i,k,j]=float(l_s[1])
                k+=1
# F_mean=np.mean(F,axis=-1)
# F_std=np.std(F,axis=-1)

rho=np.exp(-F/kB/T)/d**2
rho/=np.mean(rho[:,150:],axis=1,keepdims=True)
rho_mean=np.mean(rho,axis=-1)
rho_std=np.std(rho,axis=-1)
B22=np.sum(((1-rho)*2*np.pi*d**2)[:,:-1]*(d[:,1:]-d[:,:-1]),axis=1)
B22_mean=np.mean(B22,axis=-1)
B22_std=np.std(B22,axis=-1)

"""Plot"""
fig,axs=plt.subplots(1,2,figsize=(25,10))
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
ann=['(A)','(B)']
clrs=['r','g','b']

ax=axs[0]
ax.annotate(ann[0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.fill_between(d[i,:,0],rho_mean[i]-rho_std[i],rho_mean[i]+rho_std[i],color=clrs[i],linewidth=0,alpha=0.2)
    ax.plot(d[i,:,0],rho_mean[i],color=clrs[i],linewidth=wth,label=Mdls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
# ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_yscale('log')
ax.set_xlabel('$d_\\mathrm{com}~(\\mathrm{\\AA})$',fontdict=font)
# ax.set_ylabel('$F(d_\\mathrm{com})$ (kcal/mol) at 300 K',fontdict=font)
ax.set_ylabel('$g(d_\\mathrm{com})$ at 300 K',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ax=axs[1]
ax.annotate(ann[1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.bar(i,-np.log10(-B22_mean[i]),width=0.3,color=clrs[i],edgecolor='k',linewidth=wth)
    ax.errorbar([i],-np.log10(-B22_mean[i])[np.newaxis],yerr=[-np.log10(B22_mean[i]/(B22_mean[i]-B22_std[i]))[np.newaxis],-np.log10((B22_mean[i]+B22_std[i])/B22_mean[i])[np.newaxis]],
                        linestyle='',marker='.',color='k',markersize=5*wth,markeredgewidth=wth,
                        ecolor='k',elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nmd),Mdls,fontdict=font)
ax.set_yticks([0,-5,-10,-15],['$-10^{0}$','$-10^{5}$','$-10^{10}$','$-10^{15}$'],fontdict=font)
ax.set_ylabel('$B_{22}~(\\mathrm{\\AA^3})$',fontdict=font)

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

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()