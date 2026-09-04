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
file1='exp_cont_prob_fus.pyw'
file2='coeff_20/sim_cont_prob_fus.pyw'
figname='cont_prob_itr'+file2.split('/')[-2].split('_')[-1]
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res)
nat=len(seq)

"""Read Data 1"""
P1_tmp=pyw(file1,'Probability')[0]
P1=np.zeros((nat,nat))
for i in range(nat):
    for j in range(nat):
        P1[i,j]=P1_tmp[i*nat+j]

Pl1=np.zeros(nat)
for i in range(nat):
    for j in range(nat-i):
        Pl1[i]+=P1[i+j,j]
Pl1/=np.arange(nat,0,-1)
dsep=np.arange(nat)

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P1[i,j]=0
p1=np.zeros((ntp,ntp))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p1[res.index(seq[i]),res.index(seq[j])]+=P1[i,j]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p1[res.index(seq[j]),res.index(seq[i])]+=P1[i,j]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p1[i,j]/=num[i,j]

pl1=np.sum(p1,axis=0)
xtp=np.arange(ntp)

"""Read Data 2"""
P2_tmp=pyw(file2,'Probability')[0]
P2=np.zeros((nat,nat))
for i in range(nat):
    for j in range(nat):
        P2[i,j]=P2_tmp[i*nat+j]

Pl2=np.zeros(nat)
for i in range(nat):
    for j in range(nat-i):
        Pl2[i]+=P2[i+j,j]
Pl2/=np.arange(nat,0,-1)
dsep=np.arange(nat)

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P2[i,j]=0
p2=np.zeros((ntp,ntp))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p2[res.index(seq[i]),res.index(seq[j])]+=P2[i,j]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p2[res.index(seq[j]),res.index(seq[i])]+=P2[i,j]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p2[i,j]/=num[i,j]

pl2=np.sum(p2,axis=0)
xtp=np.arange(ntp)

"""Plot"""
plt.figure(figsize=(15,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
lenbar=8
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

# norm=colors.LogNorm(vmin=1e-4,vmax=1)
# ax=plt.subplot(2,2,1)
# img=ax.imshow(P1,cmap=cmap,norm=norm)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xlabel('Residue Index',fontsize=size)
# ax.set_ylabel('Residue Index',fontsize=size)
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# cbar=plt.colorbar(img)
# cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
# cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
# cbar.outline.set_linewidth(wth)
# cbar.set_label(f'Contact Probability',fontsize=size)

# ax=plt.subplot(2,2,2)
# norm=colors.LogNorm(vmin=np.min(p1[p1!=0]),vmax=np.max(p1))
# img=ax.imshow(p1,cmap=cmap,norm=norm)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xlabel('Residue Type',fontsize=size)
# ax.set_ylabel('Residue Type',fontsize=size)
# plt.xticks(xtp,res)
# plt.yticks(xtp,res)
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# cbar=plt.colorbar(img)
# cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
# cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
# cbar.outline.set_linewidth(wth)
# cbar.set_label(f'Contact Probability',fontsize=size)

ax=plt.subplot(1,3,1)
ax.plot(P1.reshape(-1),P2.reshape(-1),'.',linewidth=wth)
ax.plot([np.min(P2[P2!=0]),np.max(P2)],[np.min(P2[P2!=0]),np.max(P2)],linewidth=wth)
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
ax.set_xlabel('Contact Probability\nin AA Model',fontsize=size)
ax.set_ylabel('Contact Probability\nin ME Model',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')

ax=plt.subplot(1,3,2)
ax.plot(P1.reshape(-1),P2.reshape(-1),'.',linewidth=wth)
ax.plot([np.min(P2[P2!=0]),np.max(P2)],[np.min(P2[P2!=0]),np.max(P2)],linewidth=wth)
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
ax.set_xlabel('Contact Probability\nin AA Model',fontsize=size)
ax.set_ylabel('Contact Probability\nin ME Model',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')

ax=plt.subplot(1,3,3)
ax.plot(p1.reshape(-1),p2.reshape(-1),'.',linewidth=wth)
ax.plot([np.min(p1[p1!=0]),np.max(p1)],[np.min(p1[p1!=0]),np.max(p1)],linewidth=wth)
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
# ax.set_xlabel('Contact Probability Sum\nin AA Model',fontsize=size)
# ax.set_ylabel('Contact Probability Sum\nin ME Model',fontsize=size)
ax.set_xlabel('Contact Probability Mean\nin AA Model',fontsize=size)
ax.set_ylabel('Contact Probability Mean\nin ME Model',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()