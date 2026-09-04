"""Import Modules"""
import random
from scipy.optimize import curve_fit
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
file='../../chn_clst_annl/proc.npz'
figname='chn_clst'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
nrep=20
nfr=10000
nch=50

"""Read Data"""
q=np.load(file)
t=q['t']
nclst=q['nclst']
nclst_mean=q['nclst_mean']
nclst_std=q['nclst_std']
treq=q['treq']
treq_mean=q['treq_mean']
treq_std=q['treq_std']
nisol=q['nisol']
nisol_mean=q['nisol_mean']
nisol_std=q['nisol_std']
nsq=q['nsq']
dnsq=q['dnsq']
dnsq_mean=q['dnsq_mean']
dnsq_std=q['dnsq_std']
dnmin=q['dnmin']
dnmin_mean=q['dnmin_mean']
dnmin_std=q['dnmin_std']
ndnmin=np.sum(dnmin!=0,axis=-1)
dtreq=treq[:,:,1:]-treq[:,:,:-1]
dtreq_mean=np.mean(dtreq,axis=1)
# dtreq_mean=treq_mean[:,1:]-treq_mean[:,:-1]

"""Fit nclst"""
def exp1_nclst(x,a):
    return (nch-1)*np.exp(-a*x)+1
Par1_nclst=np.ones(nmd)
for i in range(nmd):
    par,cov=curve_fit(exp1_nclst,t,nclst_mean[i])
    Par1_nclst[i]=par[0]

def exp2_nclst(x,A,a,b):
    return (nch-1)*(A*np.exp(-a*x)+(1-A)*np.exp(-b*x))+1
init=[0.5,1.0,1.0]
bounds=([0.0,0.0,0.0],
        [1.0,np.inf,np.inf])
Par2_nclst=np.ones((nmd,3))
for i in range(nmd):
    par,cov=curve_fit(exp2_nclst,t,nclst_mean[i],p0=init,bounds=bounds)
    Par2_nclst[i]=par
    if Par2_nclst[i,1]<Par2_nclst[i,2]:
        Par2_nclst[i,0]=1-Par2_nclst[i,0]
        Par2_nclst[i,1:]=Par2_nclst[i,1:][::-1]

S0=nch
Sinf=1
P=Par2_nclst[:,0]
k1=Par2_nclst[:,1]
k2=Par2_nclst[:,2]
alf=S0
bet=((P*k1+(1-P)*k2)*Sinf-(1-P)*(k2-k1)*S0)/k1
gma=Sinf
print(bet)

nclst_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    nclst_fit[i]=exp2_nclst(t,Par2_nclst[i,0],Par2_nclst[i,1],Par2_nclst[i,2])
treq_fit=np.zeros((nmd,nch))
for i in range(nmd):
    ks=0
    for j in range(nch):
        for k in range(ks,nfr-1):
            if nclst_fit[i,k+1]<=nch-j<=nclst_fit[i,k]:
                treq_fit[i,nch-j-1]=(nch-j-nclst_fit[i,k])*(t[k+1]-t[k])/(nclst_fit[i,k+1]-nclst_fit[i,k])+t[k]
                ks=k
                break
dtreq_fit=treq_fit[:,1:]-treq_fit[:,:-1]

"""Fit nisol"""
def exp1_nisol(x,a):
    return nch*np.exp(-a*x)
Par1_nisol=np.ones(nmd)
for i in range(nmd):
    par,cov=curve_fit(exp1_nisol,t,nisol_mean[i])
    Par1_nisol[i]=par[0]

def exp2_nisol(x,A,a,b):
    return nch*(A*np.exp(-a*x)+(1-A)*np.exp(-b*x))
init=[0.5,1.0,1.0]
# bounds=([0.0,0.0,0.0],
#         [1.0,np.inf,np.inf])
bounds=([0.0,0.0,0.0],
        [1.0,1e1,1e1])
Par2_nisol=np.ones((nmd,3))
for i in range(nmd):
    # par,cov=curve_fit(exp2_nisol,t,nisol_mean[i],p0=init,bounds=bounds)
    # Par2_nisol[i]=par
    def loss(vars):
        A,a,b=vars
        return np.sum((exp2_nisol(t,A,a,b)-nisol_mean[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par2_nisol[i]=par
    if Par2_nisol[i,1]<Par2_nisol[i,2]:
        Par2_nisol[i,0]=1-Par2_nisol[i,0]
        Par2_nisol[i,1:]=Par2_nisol[i,1:][::-1]

S0=nch
Sinf=0
P=Par2_nisol[:,0]
k1=Par2_nisol[:,1]
k2=Par2_nisol[:,2]
alf=S0
bet=((P*k1+(1-P)*k2)*Sinf-(1-P)*(k2-k1)*S0)/k1
gma=Sinf
print(bet)

# def exp1(x,a):
#     return np.log(nch*np.exp(-a*np.exp(x)))
# Par=np.ones(nmd)
# for i in range(nmd):
#     par,cov=curve_fit(exp1,np.log(t[1:]),np.log(nclst_mean[i,1:]))
#     Par[i]=par[0]

"""Plot"""
fig,axs=plt.subplots(4,4,figsize=(25,25))
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
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nrep)
ann=[['(A)','(B)','(C)','(D)'],['(F)','(G)','(H)','(I)'],['(J)','(K)','(L)','(M)'],['(N)','(O)','(P)','(E)']]
clrs=['r','g','b']

for i in range(nmd):
    ax=axs[0,i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nrep):
        ax.plot(t,nclst[i,j],color=cmap(norm(j)),linewidth=wth,alpha=0.2)
        ax.plot(t,nclst_mean[i],color='k',linewidth=wth)
    if i==0:
        ax.plot([],[],color='k',linewidth=wth,label='Mean')
    ax.plot(t,exp1_nclst(t,Par1_nclst[i]),color=(0.7,0.2,0.0),linestyle=':',linewidth=wth,label='$y=49\\mathrm{e}^{-%.2fx}+1$'%Par1_nclst[i])
    ax.plot(t,exp2_nclst(t,Par2_nclst[i,0],Par2_nclst[i,1],Par2_nclst[i,2]),color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
            label='$y=49(%.2f\\mathrm{e}^{-%.2fx}+$\n$%.2f\\mathrm{e}^{-%.3fx})+1$'%(Par2_nclst[i,0],Par2_nclst[i,1],1-Par2_nclst[i,0],Par2_nclst[i,2]))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(-2.5,52.5)
    ax.set_xlabel('$t$ (ns) of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$N_\\mathrm{clst}$',fontdict=font)
    ax.legend(loc='upper right',prop={'family':font['family'],'size':0.5*font['size']})

ax=axs[0,3]
ax.annotate(ann[0][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t,nclst_mean[i],color=clrs[i],linewidth=wth,label=Mdls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(25))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim(-2.5,52.5)
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}N_\\mathrm{clst}\\rangle$',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ylims=[axs[0,i].get_ylim() for i in range(4)]
for i in range(4):
    axs[0,i].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nmd):
    ax=axs[1,i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nrep):
        ax.plot(np.arange(1,nch+1),treq[i,j],color=cmap(norm(j)),linewidth=wth,alpha=0.5)
    ax.plot(np.arange(1,nch+1),treq_mean[i],color='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(-2.5,52.5)
    ax.set_xlabel('$N_\\mathrm{clst}$ of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$t$ (ns)',fontdict=font)

ax=axs[1,3]
ax.annotate(ann[1][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(np.arange(1,nch+1),treq_mean[i],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_ylim(-2.5,52.5)
ax.set_xlabel('$N_\\mathrm{clst}$',fontdict=font)
ax.set_ylabel('$\\langle{}t\\rangle$ (ns)',fontdict=font)

ylims=[axs[1,i].get_ylim() for i in range(4)]
for i in range(4):
    axs[1,i].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nmd):
    ax=axs[2,i]
    ax.annotate(ann[2][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nrep):
        ax.plot(t,nisol[i,j],color=cmap(norm(j)),linewidth=wth,alpha=0.2)
    ax.plot(t,nisol_mean[i],color='k',linewidth=wth)
    ax.plot(t,exp1_nisol(t,Par1_nisol[i]),color=(0.7,0.2,0.0),linestyle=':',linewidth=wth,label='$y=50\\mathrm{e}^{-%.2fx}$'%Par1_nisol[i])
    ax.plot(t,exp2_nisol(t,Par2_nisol[i,0],Par2_nisol[i,1],Par2_nisol[i,2]),color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
            label='$y=50(%.2f\\mathrm{e}^{-%.2fx}+$\n$%.2f\\mathrm{e}^{-%.2fx})$'%(Par2_nisol[i,0],Par2_nisol[i,1],1-Par2_nisol[i,0],Par2_nisol[i,2]))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(-2.5,52.5)
    ax.set_xlabel('$t$ (ns) of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$N_\\mathrm{iso}$',fontdict=font)
    ax.legend(loc='upper right',prop={'family':font['family'],'size':0.5*font['size']})

ax=axs[2,3]
ax.annotate(ann[2][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t,nisol_mean[i],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(25))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim(-2.5,52.5)
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}N_\\mathrm{iso}\\rangle$',fontdict=font)

ylims=[axs[3,i].get_ylim() for i in range(4)]
for i in range(4):
    axs[3,i].set_ylim(np.min(ylims),np.max(ylims))

for i in range(nmd):
    ax=axs[3,i]
    ax.annotate(ann[3][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nrep):
        ax.plot(t,np.sqrt(nsq[i,j]),color=cmap(norm(j)),linewidth=wth,alpha=0.5)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(4.8,52.2)
    ax.set_xlabel('$t$ (ns) of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$A_\\mathrm{clst}$',fontdict=font)

ax=axs[3,3]
ax.annotate(ann[3][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(treq_fit[i,1:],-1/dtreq_fit[i],color=clrs[i],linewidth=wth)
ax.plot(t,[np.max(-1/dtreq_fit)]*nfr,linewidth=0)
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
ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_ylim(-2.5,52.5)
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}v_\\mathrm{clst}\\rangle$ (ns$^{-1}$)',fontdict=font)

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