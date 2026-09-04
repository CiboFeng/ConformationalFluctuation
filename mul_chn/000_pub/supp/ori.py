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
dir1='../../dst/dst'
dir2='../../ori/ori'
dir3='../../ori/ori_rel'
dir4='../../ori/mag'
figname='ori'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi_z=250
nhi_tht=50
nhi_phi=50
nhi_mag=50
maglim=[[-2.5,2.5],[-2.5,2.5],[-2.5,2.5],[0,2.5],[0,2.5],[0,2.5],[0.0,5.0]]
nq=len(maglim)

"""Read Data"""
z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]

z=np.zeros((nmd,nT,nhi_z))
tht=np.zeros((nmd,nT,nhi_tht))
ptht=np.zeros((nmd,nT,nhi_z,nhi_tht))
tht_mean=np.zeros((nmd,nT,nhi_z))
tht_std=np.zeros((nmd,nT,nhi_z))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir2}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        tht[i,j]=q['tht']
        ptht[i,j]=q['p']
        tht_mean[i,j]=q['tht_mean']
        tht_std[i,j]=q['tht_std']
tht*=180/np.pi
ptht/=180/np.pi
tht_mean*=180/np.pi
tht_std*=180/np.pi

ptht=ptht/np.sum(ptht*180/nhi_tht,axis=-1,keepdims=True)
p_rand=np.pi/360*np.sin(np.pi*tht[0,0]/180)
dptht=2*np.arccos(np.clip(np.sum(np.sqrt(ptht)*np.sqrt(p_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_tht,axis=-1),-1.0,1.0))

z=np.zeros((nmd,nT,nhi_z))
phi=np.zeros((nmd,nT,nhi_phi))
pphi=np.zeros((nmd,nT,nhi_z,nhi_phi))
phi_mean=np.zeros((nmd,nT,nhi_z))
phi_std=np.zeros((nmd,nT,nhi_z))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir3}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        phi[i,j]=q['tht']
        pphi[i,j]=q['p']
        phi_mean[i,j]=q['tht_mean']
        phi_std[i,j]=q['tht_std']
phi*=180/np.pi
pphi/=180/np.pi
phi_mean*=180/np.pi
phi_std*=180/np.pi
pphi=pphi/np.sum(pphi*180/nhi_phi,axis=-1,keepdims=True)

dpphi=2*np.arccos(np.clip(np.sum(np.sqrt(pphi)*np.sqrt(p_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_phi,axis=-1),-1.0,1.0))

z=np.zeros((nmd,nT,nhi_z))
mag=np.zeros((nmd,nT,nhi_mag,3))
pmag=np.zeros((nmd,nT,nhi_z,nhi_mag,3))
mag_mean=np.zeros((nmd,nT,nhi_z,3))
mag_std=np.zeros((nmd,nT,nhi_z,3))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir4}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        mag[i,j]=q['mag']
        pmag[i,j]=q['p']
        mag_mean[i,j]=q['mag_mean']
        mag_std[i,j]=q['mag_std']
pmag=pmag/np.sum(pmag*2/nhi_mag,axis=-2,keepdims=True)

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(25,18))
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
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-5,vmax=5e-2)
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),tht[i,j],ptht[i,j].T,cmap=cmap,norm=norm)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),tht_mean[i,j],'c',linewidth=wth,label='Mean')
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
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[j],fontdict=font)
        if j==0:
            ax.set_ylabel('$\\theta~(\\circ)$ of %s'%Mdls[i],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$f(\\theta|z_\\mathrm{norm})$',fontdict=font)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
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

fig,axs=plt.subplots(1,nT,figsize=(25,10))
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
ann=['(A)','(B)','(C)','(D)']
clrs=['r','g','b']

for i in range(nT):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),dptht[j,i],color=clrs[j],linewidth=wth,label=Mdls[j])
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
    ax.set_xlim(-1,1)
    ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$d_\\mathrm{FR}(f(\\theta),f_\\mathrm{rand}(\\theta))$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

ylims=[axs[i].get_ylim()for i in range(nT)]
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

fig,axs=plt.subplots(1,nmd,figsize=(25,10))
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
ann=['(A)','(B)','(C)']

for i in range(nmd):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),dptht[i,j],color=cmap(j/(nT-1)),linewidth=wth,label='%i K'%Ts[j])
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
    ax.set_xlim(-1,1)
    ax.set_xlabel('$z_\\mathrm{norm}$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$d_\\mathrm{FR}(f(\\theta),f_\\mathrm{rand}(\\theta))$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

ylims=[axs[i].get_ylim()for i in range(nmd)]
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

fig,axs=plt.subplots(nmd,nT,figsize=(25,18))
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
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=1)
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),phi[i,j],pphi[i,j].T*100,cmap=cmap,norm=norm)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),phi_mean[i,j],'c',linewidth=wth,label='Mean')
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
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[j],fontdict=font)
        if j==0:
            ax.set_ylabel('$\\phi~(\\circ)$ of %s'%Mdls[i],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$f(\\phi|z_\\mathrm{norm})\\times10^{-2}$',fontdict=font)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
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

plt.savefig(f'{figname}_rel.png',format='png',dpi=50)
plt.savefig(f'{figname}_rel.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nT,figsize=(25,10))
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
ann=['(A)','(B)','(C)','(D)']
clrs=['r','g','b']

for i in range(nT):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),dpphi[j,i],color=clrs[j],linewidth=wth,label=Mdls[j])
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
    ax.set_xlim(-1,1)
    ax.set_yscale('log')
    ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$d_\\mathrm{FR}(f(\\phi),f_\\mathrm{rand}(\\phi))$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

ylims=[axs[i].get_ylim()for i in range(nT)]
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

plt.savefig(f'{figname}_rel_1.png',format='png',dpi=50)
plt.savefig(f'{figname}_rel_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,nmd,figsize=(25,10))
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
ann=['(A)','(B)','(C)']

for i in range(nmd):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),dpphi[i,j],color=cmap(j/(nT-1)),linewidth=wth,label='%i K'%Ts[j])
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
    ax.set_xlim(-1,1)
    ax.set_yscale('log')
    ax.set_xlabel('$z_\\mathrm{norm}$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$d_\\mathrm{FR}(f(\\phi),f_\\mathrm{rand}(\\phi))$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

ylims=[axs[i].get_ylim()for i in range(nmd)]
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

plt.savefig(f'{figname}_rel_2.png',format='png',dpi=50)
plt.savefig(f'{figname}_rel_2.pdf',format='pdf')
# plt.show()

figlabels=['x','y','z']
ylabels=['$\\chi_x$','$\\chi_y$','$\\chi_z$']
norms=[colors.Normalize(vmin=0,vmax=1),
       colors.Normalize(vmin=0,vmax=1),
       colors.LogNorm(vmin=1e-2,vmax=2e0)]
for i in range(3):
    fig,axs=plt.subplots(nmd,nT,figsize=(25,18))
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
    ytick=1
    ctick=1
    color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
    nodes=[0/3,1/3,2/3,3/3]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=norms[i]
    ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

    for j in range(nmd):
        for k in range(nT):
            ax=axs[j,k]
            ax.annotate(ann[j][k],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
            ax.pcolormesh(z[j,k]/((z2[j,k]-z1[j,k])/2+d[j,k]),mag[j,k,:,i],pmag[j,k,:,:,i].T,cmap=cmap,norm=norm)
            ax.plot(z[j,k]/((z2[j,k]-z1[j,k])/2+d[j,k]),mag_mean[j,k,:,i],'c',linewidth=wth,label='Mean')
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
            ax.set_xlim(-1,1)
            if j==nmd-1:
                ax.set_xlabel('$z_\\mathrm{norm}$ at %d K'%Ts[k],fontdict=font)
            if k==0:
                ax.set_ylabel('%s of %s'%(ylabels[i],Mdls[j]),fontdict=font)
            if (j,k)==(0,0):
                leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
                for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                    txt.set_color(hnd.get_color())

    top=axs[0,-1].get_position().y1
    bot=axs[-1,-1].get_position().y0
    cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
    sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm,cax=cbar_ax)
    cbar.set_label('$f($%s$|z_\\mathrm{norm})$'%ylabels[i],fontdict=font)
    # cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
    # cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
    cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
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
    for j in range(nrows):
        for k in range(ncols):
            ax=axs_2d[j,k]
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
            outlfts[j,k]=boxlft-alllft
            outrigs[j,k]=allrig-boxrig
            outbots[j,k]=boxbot-allbot
            outtops[j,k]=alltop-boxtop

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

    plt.savefig(f'{figname}_mag_{figlabels[i]}.png',format='png',dpi=50)
    plt.savefig(f'{figname}_mag_{figlabels[i]}.pdf',format='pdf')
    # plt.show()

fig,axs=plt.subplots(3,nT,figsize=(25,18))
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
xtick=1
ytick=1
ctick=1
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]
clrs=['r','g','b']
ylabels=['x','y','z']

for i in range(nT):
    for j in range(3):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nmd):
            if (i,j)==(0,0):
                ax.plot([],[],color=clrs[k],linewidth=wth,label=Mdls[k])
            ax.plot(z[k,i]/((z2[k,i]-z1[k,i])/2+d[k,i]),mag_mean[k,i,:,j],color=clrs[k],linewidth=wth)
            ax.plot(z[k,i]/((z2[k,i]-z1[k,i])/2+d[k,i]),mag_std[k,i,:,j],color=clrs[k],linewidth=wth,linestyle='--')
        if (i,j)==(1,0):
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
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        if j==2:
            ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[i],fontdict=font)
        if i==0:
            ax.set_ylabel('$\\chi_{%s}$'%ylabels[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())
        if (i,j)==(1,0):
            leg=ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

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

plt.savefig(f'{figname}_mag_1.png',format='png',dpi=50)
plt.savefig(f'{figname}_mag_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(3,nmd,figsize=(25,18))
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
xtick=1
ytick=1
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)'],['(G)','(H)','(I)']]
ylabels=['x','y','z']

for i in range(nmd):
    for j in range(3):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nT):
            if (i,j)==(0,0):
                ax.plot([],[],color=cmap(k/(nT-1)),linewidth=wth,label='%i K'%Ts[k])
            ax.plot(z[i,k]/((z2[i,k]-z1[i,k])/2+d[i,k]),mag_mean[i,k,:,j],color=cmap(k/(nT-1)),linewidth=wth)
            ax.plot(z[i,k]/((z2[i,k]-z1[i,k])/2+d[i,k]),mag_std[i,k,:,j],color=cmap(k/(nT-1)),linewidth=wth,linestyle='--')
        if (i,j)==(1,0):
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
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        if j==2:
            ax.set_xlabel('$z_\\mathrm{norm}$ of %s'%Mdls[i],fontdict=font)
        if i==0:
            ax.set_ylabel('$\\chi_{%s}$'%ylabels[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())
        if (i,j)==(1,0):
            leg=ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

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

plt.savefig(f'{figname}_mag_2.png',format='png',dpi=50)
plt.savefig(f'{figname}_mag_2.pdf',format='pdf')
# plt.show()