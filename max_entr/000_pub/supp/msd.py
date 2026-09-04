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
dir='../../diff/msd'
figname='msd'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
nfr=10000
nat=163
cut=20

"""Read Data"""
t=np.zeros((nmd,nfr,nat))
d=np.zeros((nmd,nfr,nat))
for i in range(nmd):
    q=pyw(f'{dir}/{mdls[i]}_dyn.pyw','Mean-squared')
    for j in range(nfr):
        for k in range(nat):
            t[i,j,k]=q[0,j*nat+k]
            d[i,j,k]=q[2,j*nat+k]
dm=np.mean(d,axis=-1)
dmax=np.mean(d[:,200:5000],axis=1)
dmid=np.sum(d[:,:10]*np.arange(nfr-1,nfr-11,-1)[np.newaxis,:,np.newaxis],axis=1)/np.sum(np.arange(nfr-1,nfr-11,-1))

"""Fit"""
def func(x,a,b):
    return a*x**b

Par=np.ones((nmd,2))
for i in range(nmd):
    par,cov=curve_fit(func,t[i,:cut,0],dm[i,:cut])
    Par[i]=par

"""Plot"""
fig,axs=plt.subplots(1,3,figsize=(25,8))
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=-0.1
cbarwth=0.8
cbarspcwth=0.5
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
lenbar=8
xtick=0
ytick=0
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
ann=['(A)','(B)','(C)']
clrs=['r','g','b']

for i in range(nmd):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nat):
        ax.plot(t[i,:,j],d[i,:,j],color=cmap(norm(j)),linewidth=wth,alpha=0.5)
    ax.plot(t[i,:,0],dm[i,:],color='k',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='Mean')
    # ax.plot(t[i,1:cut,0],3*Par[i,0]*t[i,1:cut,0]**Par[i,1],color=(0.5,0.5,0.5),linestyle=':',linewidth=wth,label='$\\langle{}\\mathrm{\\Delta}r^2\\rangle{}\\propto{}\\mathrm{\\Delta}t^{%.2f}$'%Par[i,1])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\Delta{}t$ (ns) of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$',fontdict=font)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(40,1500)
    if i==0:
        ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
    
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$i$',fontdict=font)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar.outline.set_linewidth(wth)

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
pos=cbar.ax.get_position()
boxrig=max(boxrig,(pos.x0+pos.width)*figwth)
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
cbar.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,(margbot+outbot)/fighei,cbarwth/figwth,1.0-(margtop+outtop+margbot+outbot)/fighei])

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()