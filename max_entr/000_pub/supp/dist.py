"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
files=['../../../all_atom/dist/dist.npz',
       '../../dist/dist/1.55sgm_1.0sgm_itr20.npz',
       '../../dist/dist/1.5sgm_0.82sgm_itr20.npz',
       '../../dist/dist/1.37sgm_0.6sgm_itr20.npz']
figname='dist'
nmd=len(files)-1
nat=163
nhi=50

"""Read Data and Calculate"""
dseq=np.arange(nat)
d2=np.zeros((nmd+1,nat,nat,nhi))
p_d2=np.zeros((nmd+1,nat,nat,nhi))
d2_mean=np.zeros((nmd+1,nat,nat))
d2_std=np.zeros((nmd+1,nat,nat))
d1=np.zeros((nmd+1,nat,nhi))
p_d1=np.zeros((nmd+1,nat,nhi))
d1_mean=np.zeros((nmd+1,nat))
d1_std=np.zeros((nmd+1,nat))
for i in range(nmd+1):
    q=np.load(files[i])
    d2[i]=q['d2']
    p_d2[i]=q['p_d2']
    d2_mean[i]=q['d2_mean']
    d2_std[i]=q['d2_std']
    d1[i]=q['d1']
    p_d1[i]=q['p_d1']
    d1_mean[i]=q['d1_mean']
    d1_std[i]=q['d1_std']

"""Plot"""
fig,axs=plt.subplots(2,3,figsize=(25,12))
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=-0.1
cbarwth=0.8
cbarspcwth=1.5
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
lenbar=8
xtick=0.1
ytick=20
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
ann=[['(A)','(B)','(E)'],['(C)','(D)','(F)']]
label0=['AA','Qch','Mid','Flx']
color0=['k','r','g','b']
color_=['C1','C6']

for i in range(nmd+1):
    j,k=divmod(i,2)
    ax=axs[j,k]
    ax.annotate(ann[divmod(i,2)[0]][divmod(i,2)[1]],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for l in range(4,nat,5):
        ax.plot(d1[i,l],p_d1[i,l],color=cmap(norm(l)),linewidth=wth)
        # ax.plot([d1_mean[i,l],d1_mean[i,l]],[np.min(p_d1[i,l]),np.max(p_d1[i,l])],'--',color=rgb,linewidth=wth)
    # ax.plot([],[],'--',color='k',linewidth=wth,label='Mean')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(0,45)
    if j==1:
        ax.set_xlabel('$d_{ij} (\\mathrm{\\AA})$',fontdict=font)
    ax.set_ylabel('$f(d_{ij})~(\\mathrm{\\AA}^{-1})$ of %s'%label0[i],fontdict=font)
    # if i==0:
    #     ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

ylims=[axs[divmod(i,2)].get_ylim() for i in range(nmd+1)]
for i in range(nmd+1):
    axs[divmod(i,2)].set_ylim(np.min(ylims),np.max(ylims))

for i in range(2):
    ax=axs[i,2]
    ax.annotate(ann[i][2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax2=ax.twinx()
    for j in range(nmd+1):
        ax.plot(dseq,d1_mean[j],linewidth=4*wth,color=color_[0])
        ax.plot(dseq,d1_std[j],linewidth=4*wth,color=color_[1])
    for j in range(nmd+1):
        ax.plot(dseq,d1_mean[j],linewidth=wth,color=color0[j],label=label0[j])  
        ax.plot(dseq,d1_std[j],linewidth=wth,color=color0[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'],color=color_[0],labelcolor=color_[0])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'],color=color_[0],labelcolor=color_[0])
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(0)
    ax.spines['left'].set_color(color_[0])
    ax.yaxis.label.set_color(color_[0]) 
    ax.set_xlabel('$|i-j|$',fontdict=font)
    ax.set_ylabel('$\\langle{}d_{ij}\\rangle~(\\mathrm{\\AA})$',fontdict=font)
    if i==1:
        ax.set_xscale('log')
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'],color=color_[1],labelcolor=color_[1])
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'],color=color_[1],labelcolor=color_[1])
    # ax2.xaxis.set_major_locator(MultipleLocator(xtick))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.yaxis.set_major_locator(MultipleLocator(10))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color(color_[1])
    ax2.yaxis.label.set_color(color_[1])
    ax2.set_ylim(ax.get_ylim())
    ax2.set_ylabel('$\\mathrm{std.}(d_{ij})~(\\mathrm{\\AA})$',fontdict=font)
    if i==1:
        ax2.set_xscale('log')

top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$|i-j|$',fontdict=font)
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