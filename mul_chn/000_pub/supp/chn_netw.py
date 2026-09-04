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
dir='../../chn_netw/netw'
figname='chn_netw'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=20
nq=11

""""""
x=np.zeros((nmd,nT,nhi,nq))
p=np.zeros((nmd,nT,nhi,nq))
x_mean=np.zeros((nmd,nT,nq))
x_std=np.zeros((nmd,nT,nq))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        x[i,j]=q['x']
        p[i,j]=q['p']
        x_mean[i,j]=q['x_mean']
        x_std[i,j]=q['x_std']
idx_sort=[0,2,1,3,4,5,6,7,8,9,10]
x=x[:,:,:,idx_sort]
p=p[:,:,:,idx_sort]
x_mean=x_mean[:,:,idx_sort]
x_std=x_std[:,:,idx_sort]

idx_slct=[0,1,2,3,4,5,10]
nq=len(idx_slct)
x=x[:,:,:,idx_slct]
p=p[:,:,:,idx_slct]
x_mean=x_mean[:,:,idx_slct]
x_std=x_std[:,:,idx_slct]

p[:,:,:,3:5]/=1000

"""Plot"""
fig,axs=plt.subplots(nq,nT,figsize=(25,30))
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
xtick=[25,10,0.05,0.002,0.003,2,0.5]
ytick=[0.05,0.3,50,2,1,1,10]
ctick=1
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)'],['(M)','(N)','(O)','(P)'],['(Q)','(R)','(S)','(T)'],['(U)','(V)','(W)','(X)'],['(Y)','(Z)','(AA)','(AB)']]
clrs=['r','g','b']
qs=['$N_I^\\mathrm{deg}$','$\\langle{}N_I^\\mathrm{deg}\\rangle{}_I$',
    '$g_I^\\mathrm{bc}$','$\\langle{}g_I^\\mathrm{bc}\\rangle{}_I$','$\\mathrm{std.}(g_I^\\mathrm{bc})_I$',
    '$\\mathrm{skw.}(g_I^\\mathrm{bc})_I$','$c_\\mathrm{swn}$']

xlims=[[],[],[(),(),(),(0,0.06)],[(),(),(),(0.008,0.012)],[(),(),(),(0.003,0.01)],[],[]]
for i in range(nT):
    for j in range(nq):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nmd):
            ax.plot(x[k,i,:,j],p[k,i,:,j],color=clrs[k],linewidth=wth,label=Mdls[k])
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_major_locator(MultipleLocator(xtick[j]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(ytick[j]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel(qs[j],fontdict=font)
        if (j,i) in [(2,3),(3,3),(4,3)]:
            ax.set_xlim(xlims[j][i])
        if j==nq-1:
            ax.set_xlabel('%s at %i K'%(qs[j],Ts[i]),fontdict=font)
        if i==0:
            if j in [3,4]:
                ax.set_ylabel('$f($%s$)\\times{}10^3$'%qs[j],fontdict=font)
            else:
                ax.set_ylabel('$f($%s$)$'%qs[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

for i in range(nq):
    xlims=[axs[i,j].get_xlim() for j in range(nT)]
    ylims=[axs[i,j].get_ylim() for j in range(nT)]
    for j in range(nT):
        axs[i,j].set_xlim(np.min(xlims),np.max(xlims))
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

fig,axs=plt.subplots(nq,nmd,figsize=(25,30))
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
xtick=[25,10,0.05,0.002,0.003,2,0.5]
ytick=[0.05,0.3,50,2,1,1,10]
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)'],['(G)','(H)','(I)'],['(J)','(K)','(L)'],['(M)','(N)','(O)'],['(P)','(Q)','(R)'],['(S)','(T)','(U)']]
qs=['$N_I^\\mathrm{deg}$','$\\langle{}N_I^\\mathrm{deg}\\rangle{}_I$',
    '$g_I^\\mathrm{bc}$','$\\langle{}g_I^\\mathrm{bc}\\rangle{}_I$','$\\mathrm{std.}(g_I^\\mathrm{bc})_I$',
    '$\\mathrm{skw.}(g_I^\\mathrm{bc})_I$','$c_\\mathrm{swn}$']

xlims=[[],[],[(),(),(0,0.06)],[(),(),(0.008,0.012)],[(),(),(0.003,0.01)],[],[]]
for i in range(nmd):
    for j in range(nq):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nT):
            rgb=cmap(k/(nT-1))
            ax.plot(x[i,k,:,j],p[i,k,:,j],color=rgb,linewidth=wth,label='%i K'%Ts[k])
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_major_locator(MultipleLocator(xtick[j]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(ytick[j]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel(qs[j],fontdict=font)
        if (j,i) in [(2,2),(3,2),(4,2)]:
            ax.set_xlim(xlims[j][i])
        if j==nq-1:
            ax.set_xlabel('%s of %s'%(qs[j],Mdls[i]),fontdict=font)
        if i==0:
            if j in [3,4]:
                ax.set_ylabel('$f($%s$)\\times{}10^3$'%qs[j],fontdict=font)
            else:
                ax.set_ylabel('$f($%s$)$'%qs[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

for i in range(nq):
    xlims=[axs[i,j].get_xlim() for j in range(nmd)]
    ylims=[axs[i,j].get_ylim() for j in range(nmd)]
    for j in range(nmd):
        axs[i,j].set_xlim(np.min(xlims),np.max(xlims))
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