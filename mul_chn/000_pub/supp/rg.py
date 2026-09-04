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
dir0='../../../max_entr/rg/rg_heat'
dir1='../../dst/dst'
dir='../../rg/rg_slab'
figname='rg'
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
nhirg0=200
nhiz=250
nhirg=50
rglim=[10.0,60.0]
nsms=10

"""Read Data"""
rg0=np.zeros((nmd,nT,nhirg0))
p_rg0=np.zeros((nmd,nT,nhirg0))
rg0_mean=np.zeros((nmd,nT))
rg0_std=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Probability')
        rg0[i,j]=q[0]
        p_rg0[i,j]=q[1]
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Mean')
        rg0_mean[i,j]=q[0,0]
        q=pyw(f'{dir0}/{mdls[i]}_{Ts[j]}.pyw','Deviation')
        rg0_std[i,j]=q[0,0]
rg0_flu=np.round(np.log10(rg0_std[:,0]/rg0_mean[:,0]),1)

z=np.zeros((nmd,nT,nhiz))
rg=np.zeros((nmd,nT,nhirg))
p=np.zeros((nmd,nT,nhiz,nhirg))
rg_mean=np.zeros((nmd,nT,nhiz))
rg_std=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        rg[i,j]=q['rg']
        p[i,j]=q['p']
        rg_mean[i,j]=q['rg_mean']
        rg_std[i,j]=q['rg_std']

pz=np.sum(p,axis=-1)
pz/=np.sum(pz*(np.max(z)-np.min(z))/nhiz,axis=-1,keepdims=True)
z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
pedge=np.zeros((nmd,nT,3,nhirg))
rgedge_mean=np.zeros((nmd,nT,3))
rgedge_std=np.zeros((nmd,nT,3))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]
        idxz=[[],[],[]]
        idxz[0]+=[round(nhiz/2)-1,round(nhiz/2)]
        for k in range(nhiz):
            if z1[i,j]+d[i,j]<z[i,j,k]<z2[i,j]-d[i,j]:
                idxz[1]+=[k]
            if z1[i,j]-d[i,j]<z[i,j,k]<z1[i,j]+d[i,j] or z2[i,j]-d[i,j]<z[i,j,k]<z2[i,j]+d[i,j]:
                idxz[2]+=[k]
        for k in range(len(idxz)):
            idx=(i,j,idxz[k])
            pedge[i,j,k]=np.sum(p[idx],axis=0)
            rgedge_mean[i,j,k]=np.sum(pz[idx]*rg_mean[idx])/np.sum(pz[idx])
            rgedge_std[i,j,k]=np.sqrt(np.sum(pz[idx]*rg_std[idx]**2)/np.sum(pz[idx])+ \
                    np.sum(pz[idx]*(rg_mean[idx]-rgedge_mean[i,j,k])**2)/np.sum(pz[idx]))
pedge/=np.sum(pedge*(np.max(rg)-np.min(rg))/nhirg,axis=-1,keepdims=True)

rg_dev=np.zeros((nmd,nT,3))
rg_dev[:,:,0]=rgedge_mean[:,:,1]/np.mean(rgedge_mean[:,:,1])-1
rg_dev[:,:,1]=rgedge_std[:,:,1]/np.mean(rgedge_std[:,:,1])-1
rg_dev[:,:,2]=rgedge_std[:,:,1]/rgedge_mean[:,:,1]/np.mean(rgedge_std[:,:,1]/rgedge_mean[:,:,1])-1

rg_mean[rg_mean==0.0]=np.nan
rg_std[rg_std==0.0]=np.nan

rg_mean_sms=np.zeros((nmd,nT,nhiz))
rg_std_sms=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        for k in range(nhiz):
            idx=(i,j,list(np.arange(max(0,k-round(nsms/2)),min(nhiz,k+round(nsms/2)))))
            rg_mean_sms[i,j,k]=np.nansum(pz[idx]*rg_mean[idx])/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0))
            rg_std_sms[i,j,k]=np.sqrt(np.nansum(pz[idx]*rg_std[idx]**2)/np.nansum(pz[idx]*np.heaviside(rg_std[idx],0.0))+ \
                    np.nansum(pz[idx]*(rg_mean[idx]-rg_mean_sms[i,j,k])**2)/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0)))

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
ytick=20
ctick=1
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=np.max(p))
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)']]

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        ax.annotate(ann[i][j],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        ax.pcolormesh(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg[i,j],p[i,j].T,cmap=cmap)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_mean[i,j],'c',linewidth=wth,label='Mean')
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_major_locator(MultipleLocator(xtick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(ytick))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlim(-1,1)
        if i==nmd-1:
            ax.set_xlabel('$z_\\mathrm{norm}$ at %i K'%Ts[j],fontdict=font)
        if j==0:
            ax.set_ylabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$ of %s'%Mdls[i],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$f(z_\\mathrm{norm},R_\\mathrm{g})$',fontdict=font)
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

fig,axs=plt.subplots(6,nT+1,figsize=(25,25))
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
ytick=0.0
ctick=1
ann=[['(A)','(B)','(C)','(D)','(E)'],['(F)','(G)','(H)','(I)','(J)'],['(K)','(L)','(M)','(N)','(O)'],['(P)','(Q)','(R)','(S)','(T)'],['(U)','(V)','(W)','(X)'],['(Y)','(Z)','(AA)','(AB)']]
clrs=['r','g','b']
ylabels=['Center','Inside','Surface','Isolated']
xlabels=['$\\langle{}R_\\mathrm{g}\\rangle$','$\\mathrm{std.}(R_\\mathrm{g})$','$\\mathrm{cv.}(R_\\mathrm{g})$']

for i in range(nT):
    for j in range(3):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nmd):
            ax.plot(rg[k,i],pedge[k,i,j],color=clrs[k],linewidth=wth,label=Mdls[k])
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)
        ax.set_ylabel('$f(R_\\mathrm{g})$ at %s'%ylabels[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

ylims=[axs[j,i].get_ylim() for j in range(3) for i in range(nT)]
for i in range(nT):
    for j in range(3):
        axs[j,i].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nT):
    ax=axs[3,i]
    ax.annotate(ann[3][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(rg0[j,i],p_rg0[j,i],color=clrs[j],linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)
    ax.set_ylabel('$f(R_\\mathrm{g})$ at %s'%ylabels[j],fontdict=font)

ylims=[axs[3,i].get_ylim() for i in range(nT)]
for i in range(nT):
    axs[3,i].set_ylim(np.min(ylims),np.max(ylims))

xlims=[axs[j,i].get_xlim() for j in range(4) for i in range(nT)]
for i in range(nT):
    for j in range(4):
        axs[j,i].set_xlim(np.min(xlims),np.max(xlims))

spc=3
for i in range(3):
    ax=axs[i,nT]
    ax.annotate(ann[i][nT],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.errorbar(np.arange(nT)*(nmd+spc)+j,rgedge_mean[j,:,i],yerr=rgedge_std[j,:,i],
                        linestyle='',marker='.',color=clrs[j],markersize=5*wth,markeredgewidth=wth,
                        ecolor=clrs[j],elinewidth=wth,capsize=2*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(nT)*(nmd+spc)+nmd/2,np.array(Ts).astype(int),fontdict=font)
    ax.set_ylabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)

ax=axs[3,nT]
ax.annotate(ann[3][nT],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,rg0_mean[i,:],yerr=rg0_std[i,:],
                    linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                    ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nT)*(nmd+spc)+nmd/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)

ylims=[axs[i,-1].get_ylim() for i in range(4)]
for i in range(4):
    axs[i,-1].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nT):
    ax=axs[4,i]
    ax.annotate(ann[4][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        # ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),rg_std_sms[j,i]/rg_mean_sms[j,i],color=clrs[j],linewidth=wth)
        ax.plot(z[j,i]/((z2[j,i]-z1[j,i])/2+d[j,i]),rg_std[j,i]/rg_mean[j,i],color=clrs[j],linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(-275,275)
    ax.set_xlim(-1,1)
    ax.set_ylim(0.195,0.245)
    ax.set_xlabel('$z_\\mathrm{norm}$',fontdict=font)
    ax.set_ylabel('$\\mathrm{cv.}(R_\\mathrm{g})$',fontdict=font)

ymin=1.1*np.min(rg_dev)-0.1*np.max(rg_dev)
yticks=np.array([-0.1,0.0,0.1])
for i in range(nT):
    ax=axs[5,i]
    ax.annotate(ann[5][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    spc=3
    for j in range(nmd):
        ax.bar(np.arange(3)*(nmd+spc)+j+1,rg_dev[j,i]-ymin,width=1.0,color=clrs[j],edgecolor='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(3)*(nmd+spc)+nmd/2,xlabels,rotation=30,fontdict=font)
    ax.set_yticks(yticks-ymin,yticks,fontdict=font)
    ax.set_xlabel('%i K'%Ts[i],fontdict=font)
    ax.set_ylabel('Relative Deviation',fontdict=font)

ylims=[axs[5,i].get_ylim() for i in range(nT)]
for i in range(nT):
    axs[5,i].set_ylim(np.min(ylims),np.max(ylims))

axs[4,-1].axis('off')
axs[5,-1].axis('off')

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

fig,axs=plt.subplots(6,nmd+1,figsize=(25,25))
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
ytick=0.0
ctick=1
color=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.5,0.0)]
nodes=[0/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(J)','(K)','(L)'],['(M)','(N)','(O)','(P)'],['(Q)','(R)','(S)'],['(T)','(U)','(V)']]
ylabels=['Center','Inside','Surface','Isolated']
xlabels=['$\\langle{}R_\\mathrm{g}\\rangle$','$\\mathrm{std.}(R_\\mathrm{g})$','$\\mathrm{cv.}(R_\\mathrm{g})$']

for i in range(nmd):
    for j in range(3):
        ax=axs[j,i]
        ax.annotate(ann[j][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
        for k in range(nT):
            ax.plot(rg[i,k],pedge[i,k,j],color=cmap(k/(nT-1)),linewidth=wth,label='%i K'%Ts[k])
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)
        ax.set_ylabel('$f(R_\\mathrm{g})$ at %s'%ylabels[j],fontdict=font)
        if (i,j)==(0,0):
            leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
            for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
                txt.set_color(hnd.get_color())

ylims=[axs[j,i].get_ylim() for j in range(3) for i in range(nmd)]
for i in range(nmd):
    for j in range(3):
        axs[j,i].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nmd):
    ax=axs[3,i]
    ax.annotate(ann[3][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        ax.plot(rg0[i,j],p_rg0[i,j],color=cmap(j/(nT-1)),linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)
    ax.set_ylabel('$f(R_\\mathrm{g})$ at %s'%ylabels[j],fontdict=font)

ylims=[axs[3,i].get_ylim() for i in range(nmd)]
for i in range(nmd):
    axs[3,i].set_ylim(np.min(ylims),np.max(ylims))

xlims=[axs[j,i].get_xlim() for j in range(4) for i in range(nmd)]
for i in range(nmd):
    for j in range(4):
        axs[j,i].set_xlim(np.min(xlims),np.max(xlims))

spc=3
for i in range(3):
    ax=axs[i,nmd]
    ax.annotate(ann[i][nmd],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.errorbar(np.arange(nmd)*(nT+spc)+j,rgedge_mean[:,j,i],yerr=rgedge_std[:,j,i],
                        linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                        ecolor=rgb,elinewidth=wth,capsize=2*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(nmd)*(nT+spc)+nT/2,Mdls,fontdict=font)
    ax.set_ylabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)

ax=axs[3,nmd]
ax.annotate(ann[3][nmd],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+spc)+i,rg0_mean[:,i],yerr=rg0_std[:,i],
                    linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                    ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nmd)*(nT+spc)+nT/2,Mdls,fontdict=font)
ax.set_ylabel('$R_\\mathrm{g}~(\\mathrm{\\AA})$',fontdict=font)

ylims=[axs[i,-1].get_ylim() for i in range(4)]
for i in range(4):
    axs[i,-1].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nmd):
    ax=axs[4,i]
    ax.annotate(ann[4][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nT):
        # ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_std_sms[i,j]/rg_mean_sms[i,j],color=cmap(j/(nT-1)),linewidth=wth)
        ax.plot(z[i,j]/((z2[i,j]-z1[i,j])/2+d[i,j]),rg_std[i,j]/rg_mean[i,j],color=cmap(j/(nT-1)),linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(-275,275)
    ax.set_xlim(-1,1)
    ax.set_ylim(0.195,0.245)
    ax.set_xlabel('$z_\\mathrm{norm}$',fontdict=font)
    ax.set_ylabel('$\\mathrm{cv.}(R_\\mathrm{g})$',fontdict=font)

ymin=1.1*np.min(rg_dev)-0.1*np.max(rg_dev)
yticks=np.array([-0.1,0.0,0.1])
for i in range(nmd):
    ax=axs[5,i]
    ax.annotate(ann[5][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    spc=3
    for j in range(nT):
        ax.bar(np.arange(3)*(nT+spc)+j+1,rg_dev[i,j]-ymin,width=1.0,color=cmap(j/(nT-1)),edgecolor='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(3)*(nT+spc)+nT/2,xlabels,rotation=30,fontdict=font)
    ax.set_yticks(yticks-ymin,yticks,fontdict=font)
    ax.set_xlabel(Mdls[i],fontdict=font)
    ax.set_ylabel('Relative Deviation',fontdict=font)

ylims=[axs[5,i].get_ylim() for i in range(nmd)]
for i in range(nmd):
    axs[5,i].set_ylim(np.min(ylims),np.max(ylims))

axs[4,-1].axis('off')
axs[5,-1].axis('off')

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
