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
dir1='../../dst/dst_r'
dir2='../../ori/ori_r'
dir3='../../ori/ori_r_rel'
dir4='../../ori/mag_r'
figname='ori_r'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi_r=250
nhi_tht=50
nhi_phi=50
nhi_mag=50

"""Read Data"""
rs=np.zeros(nmd)
d=np.zeros(nmd)
for i in range(nmd):
    rs[i],d[i]=pyw(f'{dir1}/{mdls[i]}.pyw','Position')[:2,0]

r=np.zeros((nmd,nhi_r))
tht=np.zeros((nmd,nhi_tht))
ptht=np.zeros((nmd,nhi_r,nhi_tht))
tht_mean=np.zeros((nmd,nhi_r))
tht_std=np.zeros((nmd,nhi_r))
for i in range(nmd):
    q=np.load(f'{dir2}/{mdls[i]}.npz')
    r[i]=q['r']
    tht[i]=q['tht']
    ptht[i]=q['p']
    tht_mean[i]=q['tht_mean']
    tht_std[i]=q['tht_std']
tht*=180/np.pi
ptht/=180/np.pi
tht_mean*=180/np.pi
tht_std*=180/np.pi
for i in range(nmd):
    idx=r[i]>rs[i]+d[i]
    ptht[i,idx]=np.nan
ptht=ptht/np.sum(ptht*180/nhi_tht,axis=-1,keepdims=True)

p_rand=np.pi/360*np.sin(np.pi*tht[0]/180)
dptht=2*np.arccos(np.clip(np.sum(np.sqrt(ptht)*np.sqrt(p_rand[np.newaxis,np.newaxis])*180/nhi_tht,axis=-1),-1.0,1.0))

r=np.zeros((nmd,nhi_r))
mag=np.zeros((nmd,nhi_mag))
pmag=np.zeros((nmd,nhi_r,nhi_mag))
mag_mean=np.zeros((nmd,nhi_r))
mag_std=np.zeros((nmd,nhi_r))
for i in range(nmd):
    q=np.load(f'{dir4}/{mdls[i]}.npz')
    r[i]=q['r']
    mag[i]=q['mag']
    pmag[i]=q['p']
    mag_mean[i]=q['mag_mean']
    mag_std[i]=q['mag_std']
pmag=pmag/np.sum(pmag*2/nhi_mag,axis=-1,keepdims=True)

"""Plot"""
fig,axs=plt.subplots(2,nmd,figsize=(25,15))
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
xtick=1
ytick=90
ctick=1
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap1=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm1=colors.LogNorm(vmin=1e-5,vmax=5e-2)
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap2=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm2=colors.LogNorm(vmin=1e-2,vmax=2e0)
ann=[['(B)','(C)','(D)'],['(F)','(G)','(H)']]

for i in range(nmd):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.pcolormesh(r[i]/(rs[i]+d[i]),tht[i],ptht[i].T,cmap=cmap1,norm=norm1)
    ax.plot(r[i]/(rs[i]+d[i]),tht_mean[i],'c',linewidth=wth,label='Mean')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(90))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(0,1)
    if i==0:
        ax.set_ylabel('$\\theta~(\\circ)$ at 300 K',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

top=axs[0,-1].get_position().y1
bot=axs[0,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap1,norm=norm1)
sm.set_array([])
cbar1=plt.colorbar(sm,cax=cbar_ax)
cbar1.set_label('$f(\\theta|r_\\mathrm{norm})$',fontdict=font)
cbar1.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar1.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
cbar1.outline.set_linewidth(wth)

for i in range(nmd):
    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.pcolormesh(r[i]/(rs[i]+d[i]),mag[i],pmag[i].T,cmap=cmap2,norm=norm2)
    ax.plot(r[i]/(rs[i]+d[i]),mag_mean[i],'c',linewidth=wth,label='Mean')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(0,1)
    ax.set_xlabel('$r_\\mathrm{norm}$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$\\chi_r$ at 300 K',fontdict=font)

top=axs[1,-1].get_position().y1
bot=axs[1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap2,norm=norm2)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('$f(\\chi_r|r_\\mathrm{norm})$',fontdict=font)
cbar2.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar2.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
cbar2.outline.set_linewidth(wth)

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
pos1=cbar1.ax.get_position()
pos2=cbar2.ax.get_position()
boxrig=max(boxrig,(pos1.x0+pos1.width)*figwth,(pos2.x0+pos2.width)*figwth)
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
cbar1.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,
                       0.5-(margtop+outtop-margbot-outbot)/fighei/2+spchei/fighei/2,
                       cbarwth/figwth,
                       0.5-(margtop+outtop+margbot+outbot)/fighei/2-spchei/fighei/2])
cbar2.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,
                       (margbot+outbot)/fighei,
                       cbarwth/figwth,
                       0.5-(margtop+outtop+margbot+outbot)/fighei/2-spchei/fighei/2])

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

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
    ax.plot(r[i]/(rs[i]+d[i]),dptht[i],color=clrs[i],linewidth=wth,label=Mdls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(0,1)
ax.set_xlabel('$r_\\mathrm{norm}$',fontdict=font)
ax.set_ylabel('$d_\\mathrm{FR}(f(\\theta),f_\\mathrm{rand}(\\theta))$ at 300 K',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ax=axs[1]
ax.annotate(ann[1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(r[i]/(rs[i]+d[i]),mag_mean[i,:],color=clrs[i],linewidth=wth)
    ax.plot(r[i]/(rs[i]+d[i]),mag_std[i,:],color=clrs[i],linewidth=wth,linestyle='--')
ax.plot([],[],color='k',linewidth=wth,label='$\\langle\\cdot\\rangle$')
ax.plot([],[],color='k',linewidth=wth,linestyle='--',label='$\\mathrm{std.}(\\cdot)$')
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(0,1)
ax.set_ylim(-1,1)
ax.set_xlabel('$r_\\mathrm{norm}$',fontdict=font)
ax.set_ylabel('$\\chi_r$',fontdict=font)
leg=ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

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

fig=plt.figure(figsize=(25,10))
gs=fig.add_gridspec(nrows=1,ncols=2,width_ratios=[1,1],wspace=0.2)
axs=[fig.add_subplot(gs[i]) for i in range(gs.ncols)]
font={'family':'Arial','size':30}
ann=['(A)','(E)']
for i in range(2):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    if i==0:
        ax.text(0.0,0.5,'$\\theta r$',fontfamily=font['family'],fontsize=font['size'])
    ax.set_axis_off()
figwth,fighei=fig.get_size_inches()
plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)
plt.savefig(f'{figname}_blk.png',format='png')
plt.savefig(f'{figname}_blk.pdf',format='pdf')