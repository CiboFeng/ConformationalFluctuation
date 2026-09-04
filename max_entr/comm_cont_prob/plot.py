"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
files=['../../all_atom/cont_prob_mom/cont_prob_mom_1.5sgm.npy',
       '../../hps/cont_prob_mom/cont_prob_mom_1.5sgm.npy',
       'comm_cont_prob/1.55sgm_1.0sgm_itr20.pyw',
       'comm_cont_prob/1.5sgm_0.82sgm_itr20.pyw',
       'comm_cont_prob/1.37sgm_0.6sgm_itr20.pyw']
figname='comm_cont_prob'
nmd=len(files)-2
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res)
nat=len(seq)

"""Read Contact Probability"""
P=np.zeros((nmd+2,nat,nat))
P[0]=np.load(files[0])[:,:,0]
P[1]=np.load(files[1])[:,:,0]
for i in range(2,nmd+2):
    P_tmp=pyw(files[i],'Contact')
    for j in range(nat):
        for k in range(nat):
            P[i,j,k]=P_tmp[0,j*nat+k]

"""Normalize Contact Probability by Sequence Distance"""
pd=np.zeros((nmd+2,nat))
for i in range(nat):
    for j in range(nat-i):
        pd[:,i]+=P[:,j,i+j]
    pd[:,i]/=nat-i
P_=np.zeros((nmd+2,nat,nat))
for i in range(nat):
    for j in range(nat):
        P_[:,i,j]=P[:,i,j]/pd[:,abs(i-j)]

"""Average Contact Probability by Residue Type"""
for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P[:,i,j]=0
p=np.zeros((nmd+2,ntp,ntp))
num=np.zeros((nmd+2,ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p[:,res.index(seq[i]),res.index(seq[j])]+=P[:,i,j]
        num[:,res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p[:,res.index(seq[j]),res.index(seq[i])]+=P[:,i,j]
            num[:,res.index(seq[j]),res.index(seq[i])]+=1
for i in range(nmd+2):
    for j in range(ntp):
        for k in range(ntp):
            if num[i,j,k]!=0:
                p[i,j,k]/=num[i,j,k]

"""Average Normalized Contact Probability by Residue Type"""
p_=np.zeros((nmd+2,ntp,ntp))
num=np.zeros((nmd+2,ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p_[:,res.index(seq[i]),res.index(seq[j])]+=P_[:,i,j]
        num[:,res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p_[:,res.index(seq[j]),res.index(seq[i])]+=P_[:,i,j]
            num[:,res.index(seq[j]),res.index(seq[i])]+=1
for i in range(nmd+2):
    for j in range(ntp):
        for k in range(ntp):
            if num[i,j,k]!=0:
                p_[i,j,k]/=num[i,j,k]

"""Plot"""
plt.figure(figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
lenbar=8
labels=['AA','HPS','Qch','Mid','Flx']
colors=['k',(0.5,0.5,0.5),'r','g','b']

ax=plt.subplot(2,3,1)
for i in range(1,nmd+2):
    ax.plot(P[0].reshape(-1),P[i].reshape(-1),'.',color=colors[i],linewidth=wth,label=f'{labels[i]}')
ax.plot([np.min(P[0]),np.max(P[0])],[np.min(P[0]),np.max(P[0])],color='k',linewidth=wth)
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
ax.set_xlabel(f'Contact Probability\nof AA Model',fontsize=size)
ax.set_ylabel(f'Contact Probability\nof CG Model',fontsize=size)
# ax.legend(loc='best',fontsize=size)

ax=plt.subplot(2,3,2)
for i in range(1,nmd+2):
    ax.plot(P[0].reshape(-1),P[i].reshape(-1),'.',color=colors[i],linewidth=wth,label=f'{labels[i]}')
ax.plot([np.min(P[0][P[0]!=0]),np.max(P[0])],[np.min(P[0][P[0]!=0]),np.max(P[0])],color='k',linewidth=wth)
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
ax.set_xlabel(f'Contact Probability\nof AA Model',fontsize=size)
ax.set_ylabel(f'Contact Probability\nof CG Model',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='best',fontsize=size)

ax=plt.subplot(2,3,3)
for i in range(1,nmd+2):
    ax.plot(P_[0].reshape(-1),P_[i].reshape(-1),'.',color=colors[i],linewidth=wth,label=f'{labels[i]}')
ax.plot([np.min(P_[0]),np.max(P_[0])],[np.min(P_[0]),np.max(P_[0])],color='k',linewidth=wth)
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
ax.set_xlabel(f'Distance-free Contact\nof AA Model',fontsize=size)
ax.set_ylabel(f'Distance-free Contact\nof CG Model',fontsize=size)
# ax.legend(loc='best',fontsize=size)

ax=plt.subplot(2,3,4)
ax.plot([np.min(p[0]),np.max(p[0])],[np.min(p[0]),np.max(p[0])],color='k',linewidth=wth)
for i in range(1,nmd+2):
    ax.plot(p[0].reshape(-1),p[i].reshape(-1),'.',color=colors[i],linewidth=wth,label=f'{labels[i]}')
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
ax.set_xlabel(f'Type-mean Contact Probability\nof AA Model',fontsize=size)
ax.set_ylabel(f'Type-mean Contact Probability\nof CG Model',fontsize=size)
# ax.legend(loc='best',fontsize=size)

ax=plt.subplot(2,3,6)
for i in range(1,nmd+2):
    ax.plot(p_[0].reshape(-1),p_[i].reshape(-1),'.',color=colors[i],linewidth=wth,label=f'{labels[i]}')
ax.plot([np.min(p_[0]),np.max(p_[0])],[np.min(p_[0]),np.max(p_[0])],color='k',linewidth=wth)
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
ax.set_xlabel(f'Type-mean Distance-free\nContact of AA Model',fontsize=size)
ax.set_ylabel(f'Type-mean Distance-free\nContact of CG Model',fontsize=size)
# ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()