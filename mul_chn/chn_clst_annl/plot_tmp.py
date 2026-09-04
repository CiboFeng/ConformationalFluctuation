"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
mdls=['chn_clst/1.37sgm_0.6sgm',
      'chn_clst/1.5sgm_0.82sgm',
      'chn_clst/1.55sgm_1.0sgm']
figname='chn_clst'
nrep=20
nch=50
dt=10e-6*10000 ### ns

"""Read Data and Calculate"""
nclst=[]
nsgl=[]
nsq=[]
for i in range(len(mdls)):
    nclst.append([])
    nsgl.append([])
    nsq.append([])
    for j in range(nrep):
        nclst[-1].append([])
        nsgl[-1].append([])
        nsq[-1].append([])
        q=pyw(f'{mdls[i]}_{j}.pyw','Cluster')
        nfr=round(np.shape(q)[1]/nch)
        Q=np.zeros((nfr,nch))
        for k in range(np.shape(q)[1]):
                Q[int(q[0,k]),int(q[1,k])]=int(q[2,k])
        for k in range(nfr):
            nclst[-1][-1]+=[len(list(set(Q[k])))]
            num=0
            num1=0
            for l in range(int(np.max(Q[k]))+1):
                numl=np.sum(Q[k]==l)
                num+=(numl)**2
                if numl==1.0:
                    num1+=1
            nsgl[-1][-1]+=[num1]
            nsq[-1][-1]+=[num]

nfr=np.zeros((len(mdls),nrep),dtype=int)
for i in range(len(mdls)):
     for j in range(nrep):
          nfr[i,j]=len(nclst[i][j])
nfrmax=np.max(nfr)
t=np.arange(nfrmax)*dt

nclst_mean=np.zeros((len(mdls),nfrmax))
for i in range(len(mdls)):
    for j in range(nfrmax):
        nsuc=0
        for k in range(nrep):
            try:
                nclst_mean[i,j]+=nclst[i][k][j]
                nsuc+=1
            except Exception as err:
                pass
        nclst_mean[i,j]/=nsuc

"""Fit"""
def func(x,a):
    return nch*np.exp(-a*x)

Par=np.ones(len(mdls))
for i in range(len(mdls)):
    par,cov=curve_fit(func,t,nclst_mean[i])
    Par[i]=par[0]

# def func(x,a):
#     return np.log(nch*np.exp(-a*np.exp(x)))

# Par=np.ones(len(mdls))
# for i in range(len(mdls)):
#     par,cov=curve_fit(func,np.log(t[1:]),np.log(nclst_mean[i,1:]))
#     Par[i]=par[0]

"""Plot"""
fig=plt.figure(figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nrep)
clrs=['r','g','b']
lbls=['Flex','Mid','Fix']

for i in range(len(mdls)):
    ax=plt.subplot(2,2,i+1)
    for j in range(nrep):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(np.arange(len(nclst[i][j]))*dt,nclst[i][j],color=rgb,linewidth=wth,alpha=0.5)
    ax.plot(t,nclst_mean[i],color='k',linewidth=wth,label='Mean')
    ax.plot(t,nch*np.exp(-Par[i]*t),color='k',linestyle=':',linewidth=wth,label='$y=50\\mathrm{e}^{-%.2fx}$'%Par[i])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(-2.5,52.5)
    if divmod(i,2)[0]==1:
        ax.set_xlabel('Time (ns)',fontsize=size)
    ax.set_ylabel('Number of Clusters\n of %s Model'%lbls[i],fontsize=size)
    ax.legend(loc='upper right',fontsize=size)

ax=plt.subplot(2,2,4)
for i in range(len(mdls)):
    ax.plot(t,nclst_mean[i],color=clrs[i],linewidth=wth,label=lbls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim(-2.5,52.5)
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Number of Clusters',fontsize=size)
ax.legend(loc='upper right',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_num.png',format='png')
plt.savefig(f'{figname}_num.pdf',format='pdf')
plt.show()

fig=plt.figure(figsize=(15,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nrep)
lbls=['Flex','Mid','Fix']

for i in range(len(mdls)):
    ax=plt.subplot(1,3,i+1)
    for j in range(nrep):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(np.arange(len(nsgl[i][j]))*dt,nsgl[i][j],color=rgb,linewidth=wth,alpha=0.5)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_ylim(-2.5,52.5)
    ax.set_xlabel('Time (ns) of %s Model'%lbls[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Number of Sparate Chain',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_sgl.png',format='png')
plt.savefig(f'{figname}_sgl.pdf',format='pdf')
plt.show()

fig=plt.figure(figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nrep)
lbls=['Flex','Mid','Fix']

for i in range(len(mdls)):
    ax=plt.subplot(2,3,i+1)
    for j in range(nrep):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(np.arange(len(nsq[i][j]))*dt,np.sqrt(nsq[i][j]),color=rgb,linewidth=wth,alpha=0.5)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(4.8,52.2)
    # ax.set_xlabel('Time (ns) of %s Model'%lbls[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Clusters Amount',fontsize=size)

for i in range(len(mdls)):
    ax=plt.subplot(2,3,len(mdls)+i+1)
    for j in range(nrep):
        normc=norm(j)
        rgb=cmap(normc)
        nsq_=np.array(nsq[i][j])
        dn=np.sqrt((nsq_[1:]-nsq_[:-1])/2)
        dn[dn==0.0]=np.nan
        ax.scatter(np.arange(len(dn))*dt,dn,color=rgb,s=10*wth,alpha=0.5)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_ylim(-1.6,32.6)
    ax.set_xlabel('Time (ns) of %s Model'%lbls[i],fontsize=size)
    if i==0:
        ax.set_ylabel('Merging Amount',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_amt.png',format='png')
plt.savefig(f'{figname}_amt.pdf',format='pdf')
plt.show()