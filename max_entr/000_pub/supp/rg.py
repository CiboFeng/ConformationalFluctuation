"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
file0='../../../all_atom/rg/rg.pyw'
file='../../rg/fus_sgm/rg.pyw'
figname='rg'
Mdls=['Qch','Mid','Flx']
nhi=200
csgm=[[1.5,1.55,1.6],
      [1.5,1.5,1.5,1.5],
      [1.37,1.4,1.5]]
esgm=[[1.0,1.0,1.0],
      [0.6,0.8,0.82,1.0],
      [0.6,0.6,0.6]]
itr=[[20,20,20],
     [20,20,20,20],
     [20,20,20]]
par=[[1.5,1.55,1.6],
     [0.6,0.8,0.82,1.0],
     [1.37,1.4,1.5]]
par0=[1.0,1.5,0.6]

"""Read Data and Calculate"""
q=pyw(file0,'Probability')
rg0=q[0]
p_rg0=q[1]
q=pyw(file0,'Mean')
rg0_mean=q[0][0]
q=pyw(file0,'Deviation')
rg0_std=q[0][0]

rg=[]
p_rg=[]
rg_mean=[]
rg_std=[]
for i in range(len(par)):
    rg.append([])
    p_rg.append([])
    rg_mean.append([])
    rg_std.append([])
    for j in range(len(par[i])):
        q=pyw(file.replace('sgm',f'{csgm[i][j]}sgm_{esgm[i][j]}sgm').replace('rg.pyw',f'rg_itr{itr[i][j]}.pyw'),'Probability')
        rg[i]+=[q[0]]
        p_rg[i]+=[q[1]]
        q=pyw(file.replace('sgm',f'{csgm[i][j]}sgm_{esgm[i][j]}sgm').replace('rg.pyw',f'rg_itr{itr[i][j]}.pyw'),'Mean')
        rg_mean[i]+=[q[0][0]]
        q=pyw(file.replace('sgm',f'{csgm[i][j]}sgm_{esgm[i][j]}sgm').replace('rg.pyw',f'rg_itr{itr[i][j]}.pyw'),'Deviation')
        rg_std[i]+=[q[0][0]]

"""Plot"""
fig,axs=plt.subplots(2,3,figsize=(25,15))
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
xtick=[[],[0.05,0.2,0.05]]
ytick=0
ctick=1
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)']]
clrs=[['C0','C1','C2'],['C3','C4','C5','C6'],['C7','C8','C9']]
var=['$c_\\mathrm{c}$','$c_\\mathrm{e}$','$c_\\mathrm{c}$']
var0=['$c_\\mathrm{e}$','$c_\\mathrm{c}$','$c_\\mathrm{e}$']

for i in range(len(par)):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(len(par[i])):
        ax.plot(rg[i][j],p_rg[i][j],color=clrs[i][j],linewidth=wth,label='%s=%s'%(var[i],par[i][j]))
    ax.plot(rg0,p_rg0,'k',linewidth=wth,label='AA')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(11,21)
    ax.set_xlabel('$R_\\mathrm{g} (\\mathrm{\\AA})$ of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$f(R_\\mathrm{g}) (\\mathrm{\\AA}^{-1})$',fontdict=font)
    leg=ax.legend(loc='upper right',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
    for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
        txt.set_color(hnd.get_color())

    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.plot(par[i],rg_mean[i],'--',color=(0.5,0.5,0.5),linewidth=wth)
    ax.fill_between([par[i][0],par[i][-1]],[rg0_mean-rg0_std,rg0_mean-rg0_std],[rg0_mean+rg0_std,rg0_mean+rg0_std],color='k',linewidth=0,alpha=0.2)
    ax.plot([par[i][0],par[i][-1]],[rg0_mean,rg0_mean],'k',linewidth=wth)
    for j in range(len(par[i])):
        ax.errorbar(par[i][j],rg_mean[i][j],yerr=rg_std[i][j],
                        marker='.',markersize=0,markeredgewidth=2*wth,
                        ecolor=(0.5,0.5,0.5),elinewidth=2*wth,capsize=2.5*wth)
        ax.errorbar(par[i][j],rg_mean[i][j],yerr=rg_std[i][j],
                        marker='.',markersize=0,markeredgewidth=wth,
                        ecolor=clrs[i][j],elinewidth=wth,capsize=2*wth)
        ax.plot(par[i][j],rg_mean[i][j],'.',color=(0.5,0.5,0.5),markersize=8*wth)
        ax.plot(par[i][j],rg_mean[i][j],'.',color=clrs[i][j],markersize=5*wth)
    ax.plot([],[],label='%s=%s'%(var0[i],par0[i]))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(xtick[1][i]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel(var[i]+' of %s'%Mdls[i],fontdict=font)
    if i==0:
        ax.set_ylabel('$R_\\mathrm{g} (\\mathrm{\\AA})$',fontdict=font)
    ax.legend(loc='upper center',handlelength=0.0,prop={'family':font['family'],'size':font['size']})

for i in range(2):
    ylims=[axs[i,j].get_ylim() for j in range(len(par))]
    for j in range(len(par)):
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

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()