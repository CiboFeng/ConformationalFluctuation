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
file1='../../all_atom/cont_prob_mom/cont_prob_mom_6.0.npy'
file2='cont_prob_mom_6.0.npy'
figname='cont_prob_cor_all_atom_6.0'
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res)
nat=len(seq)

"""Read Data 1"""
P1=np.load(file1)
nmo=np.shape(P1)[-1]

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P1[i,j,:]=0
p1=np.zeros((ntp,ntp,nmo))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p1[res.index(seq[i]),res.index(seq[j]),:]+=P1[i,j,:]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p1[res.index(seq[j]),res.index(seq[i]),:]+=P1[i,j,:]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p1[i,j]/=num[i,j]

pl1=np.sum(p1,axis=0)

"""Read Data 2"""
P2=np.load(file2)
nmo=np.shape(P2)[-1]

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P2[i,j,:]=0
p2=np.zeros((ntp,ntp,nmo))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p2[res.index(seq[i]),res.index(seq[j]),:]+=P2[i,j,:]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p2[res.index(seq[j]),res.index(seq[i]),:]+=P2[i,j,:]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p2[i,j]/=num[i,j]

pl2=np.sum(p2,axis=0)

"""Plot"""
plt.figure(figsize=(8,4))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
lenbar=8

ax=plt.subplot(1,2,1)
# ax.scatter(p1[:,:,0].reshape(-1),p2[:,:,0].reshape(-1),s=wth**2)
ax.plot(p1[:,:,0].reshape(-1),p2[:,:,0].reshape(-1),'.',linewidth=wth)
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
ax.set_xlabel('All-atom',fontsize=size)
ax.set_ylabel('HPS',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')

ax=plt.subplot(1,2,2)
# ax.scatter(pl1[:,0],pl2[:,0],s=wth**2)
ax.plot(pl1[:,0],pl2[:,0],'.',linewidth=wth)
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
ax.set_xlabel('All-atom',fontsize=size)
ax.set_ylabel('HPS',fontsize=size)
xspan=np.max(pl1[:,0])-np.min(pl1[:,0][pl1[:,0]!=0])
yspan=np.max(pl2[:,0])-np.min(pl2[:,0][pl2[:,0]!=0])
ax.set_xlim(np.min(pl1[:,0][pl1[:,0]!=0])-0.1*xspan,np.max(pl1[:,0])+0.1*xspan)
ax.set_ylim(np.min(pl2[:,0][pl2[:,0]!=0])-0.1*yspan,np.max(pl2[:,0])+0.1*yspan)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()