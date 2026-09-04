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
file1='../../chn_clst_annl/proc.npz'
file2='../../cont/cont_frac_kin_annl.pyw'
figname='kin'
Mdls=['Qch','Mid','Flx']
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nrep=20
nfr=10000
nch=50
dt=10e-6*10000 ### ns

"""Read Data"""
q=np.load(file1)
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

"""Contact"""
Qintra=np.zeros((nmd,nfr))
Qinter=np.zeros((nmd,nfr))
for i in range(nmd):
    q=pyw(file2,f'({Mdls[i]})')
    Qintra[i]=q[1]
    Qinter[i]=q[2]
Qintras=Qintra[:,1:2]
Qintrae=np.mean(Qintra[:,-1000:],axis=-1,keepdims=True)
Qintra=(Qintra-Qintras)/(Qintrae-Qintras)
Qinters=Qinter[:,1:2]
Qintere=np.mean(Qinter[:,-1000:],axis=-1,keepdims=True)
Qinter=(Qinter-Qinters)/(Qintere-Qinters)

"""Fit Qintra"""
def exp(x,pars):
    y=np.sum(pars[::2])
    n=round(len(pars)/2)
    for i in range(0,2*n,2):
        y+=-pars[i]*np.exp(-pars[i+1]*x)
    return y

nexp_intra=3

init=[1.0,1.0]*nexp_intra
bounds=([-1e1,0.0]*nexp_intra,
        [1e3,1e1]*nexp_intra)
Par_intra=np.ones((nmd,2*nexp_intra))
for i in range(nmd):
    def loss(vars):
        return np.sum((exp(t,vars)-Qintra[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par_intra[i]=par
    idx=np.argsort(Par_intra[i,1::2])[::-1]
    Par_intra[i,::2]=Par_intra[i,::2][idx]
    Par_intra[i,1::2]=Par_intra[i,1::2][idx]

Qintra_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    Qintra_fit[i]=exp(t,Par_intra[i])
dQintra_fit=Qintra_fit[:,1:]-Qintra_fit[:,:-1]

"""Fit Qinter"""
nexp_inter=2

init=[1.0,1.0]*nexp_inter
bounds=([-1e1,0.0]*nexp_inter,
        [1e3,1e1]*nexp_inter)
Par_inter=np.ones((nmd,2*nexp_inter))
for i in range(nmd):
    def loss(vars):
        return np.sum((exp(t,vars)-Qinter[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par_inter[i]=par
    idx=np.argsort(Par_inter[i,1::2])[::-1]
    Par_inter[i,::2]=Par_inter[i,::2][idx]
    Par_inter[i,1::2]=Par_inter[i,1::2][idx]

Qinter_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    Qinter_fit[i]=exp(t,Par_inter[i])
dQinter_fit=Qinter_fit[:,1:]-Qinter_fit[:,:-1]

"""Plot"""
fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3)
Gs=[gs[0].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3),
    gs[1].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
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
xtick=100
ytick=[[0,1,50],[]]
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nrep)
clrs=['r','g','b']
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)']]

ax=axs[0][0]
ax.annotate(ann[0][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
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
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}N_\\mathrm{clst}\\rangle$',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

idx=[24,0]
for i in range(2):
    ax=axs[0][i+1]
    ax.annotate(ann[0][i+1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        viol=ax.violinplot([treq[j,:,idx[i]]],positions=[j],showextrema=False)
        for k,body in enumerate(viol['bodies']):
            body.set_facecolor(clrs[j])
            # body.set_edgecolor('k')
            # body.set_linewidth(wth)
            # body.set_alpha(1.0)
        ax.errorbar(j,treq_mean[j,idx[i]],yerr=treq_std[j,idx[i]],
                    linestyle='',marker='.',color=clrs[j],markersize=5*wth,markeredgewidth=wth,
                    ecolor=clrs[j],elinewidth=wth,capsize=2*wth)
        
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[0][i+1]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_yscale('log')
    ax.set_xticks(np.arange(nmd),Mdls,fontdict=font)
    ax.set_ylabel('$t_\\mathrm{req}$ (ns)\nfor %i Clusters'%(idx[i]+1),fontdict=font)

ax=axs[1][0]
ax.annotate(ann[1][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nrep):
    ax.plot(t,np.sqrt(nsq[0,i]),color=cmap(norm(i)),linewidth=wth,alpha=0.5)
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
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$A\\mathrm{_{clst}~of~FUS^{%s}}$'%Mdls[0],fontdict=font)

ax=axs[1][1]
ax.annotate(ann[1][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t,Qintra[i],color=clrs[i],linestyle='--',linewidth=wth)
    ax.plot(t,Qinter[i],color=clrs[i],linewidth=wth)
ax.plot([],[],color='k',linestyle='--',linewidth=wth,label='intra')
ax.plot([],[],color='k',linewidth=wth,label='inter')
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xscale('log')
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}Q_\\mathrm{intra}\\rangle~&~\\langle{}Q_\\mathrm{inter}\\rangle$',fontdict=font)
ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

ax=axs[1][2]
ax.annotate(ann[1][2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t[1:1001],dQintra_fit[i,:1000],color=clrs[i],linestyle='--',linewidth=wth)
    ax.plot(t[1:1001],dQinter_fit[i,:1000],color=clrs[i],linewidth=wth)
    ax.plot(t[1:],[np.max(dQinter_fit[i,:1000])]*(nfr-1),linewidth=0)
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
ax.set_xlabel('$t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle{}v_\\mathrm{intra}\\rangle~&~\\langle{}v_\\mathrm{inter}\\rangle~(\\mathrm{ns}^{-1})$',fontdict=font)

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

plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')