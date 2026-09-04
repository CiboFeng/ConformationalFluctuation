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

""""""
file='../../chn_netw/cont_chn_distr.npz'
figname='cont_chn_distr'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=100

""""""
q=np.load(file)
x=q['x']
p=q['p']

"""Plot"""
fig,axs=plt.subplots(1,nT,figsize=(25,8))
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
ytick=0.1
ctick=1
clrs=['r','g','b']
ann=['(A)','(B)','(C)','(D)']

for i in range(nT):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(x[j,i],p[j,i],color=clrs[j],linewidth=wth,label=Mdls[j])
        ax.fill_between(x[j,i],p[j,i],0.0*x[j,i],color='g',alpha=0.1)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$p_{IJ}$ at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$f(p_{IJ})$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xticks([-10.0,0.0],['$10^{-10}$','$10^{0}$'],fontdict=font)

ylims=[axs[i].get_ylim() for i in range(nT)]
for i in range(nT):
    axs[i].set_ylim(np.min(ylims),np.max(ylims))

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

fig,axs=plt.subplots(1,nmd,figsize=(25,8))
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
ytick=0.1
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=['(A)','(B)','(C)']

for i in range(nmd):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(x[i,j],p[i,j],color=cmap(j/(nT-1)),linewidth=wth,label='%i K'%Ts[j])
        ax.fill_between(x[i,j],p[i,j],0.0*x[i,j],color='g',alpha=0.1)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$p_{IJ}$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$f(p_{IJ})$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xticks([-10.0,0.0],['$10^{-10}$','$10^{0}$'],fontdict=font)

ylims=[axs[i].get_ylim() for i in range(nmd)]
for i in range(nmd):
    axs[i].set_ylim(np.min(ylims),np.max(ylims))

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