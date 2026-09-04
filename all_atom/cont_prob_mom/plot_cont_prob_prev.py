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
file_prev='prev/cont_prob_mom_1.5sgm.npy'
file='cont_prob_mom_1.5sgm.npy'
figname='cont_prob_1.5sgm_prev'
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res)
nat=len(seq)

"""Read Previous Data"""
P_prev=np.load(file_prev)
nmo=np.shape(P_prev)[-1]

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P_prev[i,j,:]=0
p_prev=np.zeros((ntp,ntp,nmo))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p_prev[res.index(seq[i]),res.index(seq[j]),:]+=P_prev[i,j,:]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p_prev[res.index(seq[j]),res.index(seq[i]),:]+=P_prev[i,j,:]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p_prev[i,j]/=num[i,j]

"""Read Data"""
P=np.load(file)
nmo=np.shape(P)[-1]

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P[i,j,:]=0
p=np.zeros((ntp,ntp,nmo))
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p[res.index(seq[i]),res.index(seq[j]),:]+=P[i,j,:]
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            p[res.index(seq[j]),res.index(seq[i]),:]+=P[i,j,:]
            num[res.index(seq[j]),res.index(seq[i])]+=1
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p[i,j]/=num[i,j]

xtp=np.arange(ntp)

"""Plot"""
plt.figure(figsize=(30,10))
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

ax=plt.subplot(2,4,1)
norm=colors.LogNorm(vmin=1e-4,vmax=1)
img=ax.imshow(P_prev[:,:,0],cmap=cmap,norm=norm)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('Residue Index',fontsize=size)
ax.set_ylabel('Residue Index',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')
cbar=plt.colorbar(img)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
cbar.set_label(f'Previous Contact Probability',fontsize=size)

ax=plt.subplot(2,4,2)
norm=colors.LogNorm(vmin=1e-4,vmax=1)
img=ax.imshow(P[:,:,0],cmap=cmap,norm=norm)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('Residue Index',fontsize=size)
ax.set_ylabel('Residue Index',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')
cbar=plt.colorbar(img)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
cbar.set_label(f'New Contact Probability',fontsize=size)

ax=plt.subplot(2,4,3)
ax.scatter(P_prev[:,:,0].reshape(-1),P[:,:,0].reshape(-1),linewidth=wth)
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
ax.set_xlabel('Previous Contact Probability',fontsize=size)
ax.set_ylabel('New Contact Probability',fontsize=size)

ax=plt.subplot(2,4,4)
ax.scatter(P_prev[:,:,0].reshape(-1),P[:,:,0].reshape(-1),linewidth=wth)
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
ax.set_xlabel('Previous Contact Probability',fontsize=size)
ax.set_ylabel('New Contact Probability',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')

ax=plt.subplot(2,4,5)
norm=colors.LogNorm(vmin=np.min(p[p!=0]),vmax=np.max(p))
img=ax.imshow(p_prev[:,:,0],cmap=cmap,norm=norm)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('Residue Type',fontsize=size)
ax.set_ylabel('Residue Type',fontsize=size)
plt.xticks(xtp,res)
plt.yticks(xtp,res)
cbar=plt.colorbar(img)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
cbar.set_label(f'Previous Mean\nContact Probability',fontsize=size)

ax=plt.subplot(2,4,6)
norm=colors.LogNorm(vmin=np.min(p[p!=0]),vmax=np.max(p))
img=ax.imshow(p[:,:,0],cmap=cmap,norm=norm)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('Residue Type',fontsize=size)
ax.set_ylabel('Residue Type',fontsize=size)
plt.xticks(xtp,res)
plt.yticks(xtp,res)
cbar=plt.colorbar(img)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
cbar.set_label(f'New Mean\nContact Probability',fontsize=size)

ax=plt.subplot(2,4,7)
ax.scatter(p_prev[:,:,0].reshape(-1),p[:,:,0].reshape(-1),linewidth=wth)
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
ax.set_xlabel('Previous Mean\nContact Probability',fontsize=size)
ax.set_ylabel('New Mean\nContact Probability',fontsize=size)

ax=plt.subplot(2,4,8)
ax.scatter(p_prev[:,:,0].reshape(-1),p[:,:,0].reshape(-1),linewidth=wth)
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
ax.set_xlabel('Previous Mean\nContact Probability',fontsize=size)
ax.set_ylabel('New Mean\nContact Probability',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()