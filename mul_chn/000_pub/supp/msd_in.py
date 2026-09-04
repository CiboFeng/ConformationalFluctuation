"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
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
file2='../../diff/proc_msd_res_in.npz'
figname='msd_in'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
nat=163
cut=[10,100]
dt=10e-6*10000 ### ns

"""Read Data"""
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

q=np.load(file2)
drm=q['dm']
drmm=q['dmm']

drmxy=np.sum(drm[:,:,:,:,:2],axis=-1)
drmz=drm[:,:,:,:,-1]
drmmxy=np.sum(drmm[:,:,:,:2],axis=-1)
drmmz=drmm[:,:,:,-1]

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
        for k in range(nat):
            lastnan(drmxy[i,j,:,k],100)
            lastnan(drmz[i,j,:,k],100)     
        lastnan(drmmxy[i,j],100)
        lastnan(drmmz[i,j],100)

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(25,18))
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
ytick=0
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nch)
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nch):
            rgb=cmap(norm(k))
            ax.plot(t,dxy[i,j,:,k],color=rgb,linewidth=wth,alpha=0.5)
            ax.plot(t,dz[i,j,:,k],color=rgb,linestyle='--',linewidth=wth,alpha=0.5)
        ax.plot(t,dmxy[i,j,:],color='k',linewidth=wth,label='$xy$')
        ax.plot(t,dmz[i,j,:],color='k',linestyle='--',linewidth=wth,label='$z$')
        ax.plot(t,[dmz[i,j,100]]*nfr,linewidth=0.0)
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
        if i==nmd-1:
            ax.set_xlabel('$\\Delta{}t$ (ns) at %s K'%int(Ts[j]),fontdict=font)
        if j==0:
            ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$ of %s'%Mdls[i],fontdict=font)
        if (i,j)==(0,0):
            ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
        ax.set_xscale('log')
        ax.set_yscale('log')

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

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nT,figsize=(25,15))
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
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)']]
clrs=['r','g','b']

for i in range(nT):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(t,dmxy[j,i],color=clrs[j],linewidth=wth,label=Mdls[j])
        ax.plot(t,dmz[j,i],color=clrs[j],linestyle='--',linewidth=wth)
    ax.plot(t,[dmz[0,i,100]]*nfr,linewidth=0.0)
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
    if i==0:
        ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(40,1500)

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(t,dmxy[j,i]/t/6,color=clrs[j],linewidth=wth)
        ax.plot(t,dmz[j,i]/t/6,color=clrs[j],linestyle='--',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='$xy$')
    ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='$z$')
    ax.plot(t,[(dmz[0,i]/t/6)[100]]*nfr,linewidth=0.0)
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
    ax.set_xlabel('$\\Delta{}t$ (ns) at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$D~(\\mathrm{\\AA}^2$/ns)',fontdict=font)
        ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
    ax.set_xscale('log')
    ax.set_yscale('log')

for i in range(2):
    ylims=[axs[i,j].get_ylim() for j in range(nT)]
    for j in range(nT):
        axs[i,j].set_ylim(np.min(ylims),np.max(ylims))

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

fig,axs=plt.subplots(2,nmd,figsize=(25,15))
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
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)']]

for i in range(nmd):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(t,dmxy[i,j],color=rgb,linewidth=wth,label='%i K'%Ts[j])
        ax.plot(t,dmz[i,j],color=rgb,linestyle='--',linewidth=wth)
    ax.plot(t,[dmz[i,0,100]]*nfr,linewidth=0.0)
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
    if i==0:
        ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(40,1500)

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(t,dmxy[i,j]/t/6,color=rgb,linewidth=wth)
        ax.plot(t,dmz[i,j]/t/6,color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='$xy$')
    ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='$z$')
    ax.plot(t,[(dmz[i,0]/t/6)[100]]*nfr,linewidth=0.0)
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
    ax.set_xlabel('$\\Delta{}t$ (ns) of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$D~(\\mathrm{\\AA}^2$/ns)',fontdict=font)
        ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
    ax.set_xscale('log')
    ax.set_yscale('log')

for i in range(2):
    ylims=[axs[i,j].get_ylim() for j in range(nmd)]
    for j in range(nmd):
        axs[i,j].set_ylim(np.min(ylims),np.max(ylims))

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
xtick=0
ytick=0
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nat):
            rgb=cmap(norm(k))
            ax.plot(t,drmxy[i,j,:,k],color=rgb,linewidth=wth,alpha=0.5)
            ax.plot(t,drmz[i,j,:,k],color=rgb,linestyle='--',linewidth=wth,alpha=0.5)
        ax.plot(t,drmmxy[i,j,:],color='k',linewidth=wth,label='$xy$')
        ax.plot(t,drmmz[i,j,:],color='k',linestyle='--',linewidth=wth,label='$z$')
        ax.plot(t,[drmmz[i,j,100]]*nfr,linewidth=0.0)
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
        if i==nmd-1:
            ax.set_xlabel('$\\Delta{}t$ (ns) at %s K'%int(Ts[j]),fontdict=font)
        if j==0:
            ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$ of %s'%Mdls[i],fontdict=font)
        if (i,j)==(0,0):
            ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
        ax.set_xscale('log')
        ax.set_yscale('log')

ylim=[np.inf,-np.inf]
for ax in axs.flat:
    y1,y2=ax.get_ylim()
    ylim[0]=min(ylim[0],y1)
    ylim[1]=max(ylim[1],y2)
for ax in axs.flat:
    ax.set_ylim(ylim[0],ylim[1])
    
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$i$',fontdict=font)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
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

plt.savefig(f'{figname}_res.png',format='png',dpi=50)
plt.savefig(f'{figname}_res.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nT,figsize=(25,15))
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
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)']]
clrs=['r','g','b']

for i in range(nT):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(t,drmmxy[j,i],color=clrs[j],linewidth=wth,label=Mdls[j])
        ax.plot(t,drmmz[j,i],color=clrs[j],linestyle='--',linewidth=wth)
    ax.plot(t,[drmmz[0,i,100]]*nfr,linewidth=0.0)
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
    if i==0:
        ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(40,1500)

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(t,drmmxy[j,i]/t/6,color=clrs[j],linewidth=wth)
        ax.plot(t,drmmz[j,i]/t/6,color=clrs[j],linestyle='--',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='$xy$')
    ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='$z$')
    ax.plot(t,[(drmmz[0,i]/t/6)[100]]*nfr,linewidth=0.0)
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
    ax.set_xlabel('$\\Delta{}t$ (ns) at %i K'%Ts[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$D~(\\mathrm{\\AA}^2$/ns)',fontdict=font)
        ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
    ax.set_xscale('log')
    ax.set_yscale('log')

for i in range(2):
    ylims=[axs[i,j].get_ylim() for j in range(nT)]
    for j in range(nT):
        axs[i,j].set_ylim(np.min(ylims),np.max(ylims))
        
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

plt.savefig(f'{figname}_res_1.png',format='png',dpi=50)
plt.savefig(f'{figname}_res_1.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(2,nmd,figsize=(25,15))
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
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)']]

for i in range(nmd):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(t,drmmxy[i,j],color=rgb,linewidth=wth,label='%i K'%Ts[j])
        ax.plot(t,drmmz[i,j],color=rgb,linestyle='--',linewidth=wth)
    ax.plot(t,[drmmz[i,0,100]]*nfr,linewidth=0.0)
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
    if i==0:
        ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA}^2)$',fontdict=font)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(40,1500)

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(t,drmmxy[i,j]/t/6,color=rgb,linewidth=wth)
        ax.plot(t,drmmz[i,j]/t/6,color=rgb,linestyle='--',linewidth=wth)
    ax.plot([],[],color='k',linewidth=wth,label='$xy$')
    ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='$z$')
    ax.plot(t,[(drmmz[i,0]/t/6)[100]]*nfr,linewidth=0.0)
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
    ax.set_xlabel('$\\Delta{}t$ (ns) of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$D~(\\mathrm{\\AA}^2$/ns)',fontdict=font)
        ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})
    ax.set_xscale('log')
    ax.set_yscale('log')

for i in range(2):
    ylims=[axs[i,j].get_ylim() for j in range(nmd)]
    for j in range(nmd):
        axs[i,j].set_ylim(np.min(ylims),np.max(ylims))

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

plt.savefig(f'{figname}_res_2.png',format='png',dpi=50)
plt.savefig(f'{figname}_res_2.pdf',format='pdf')
# plt.show()