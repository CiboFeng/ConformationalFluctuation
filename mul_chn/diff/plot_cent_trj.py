"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
file='cent_trj'
figname=file
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100

"""Read Data"""
r1=np.zeros((nmd,nT,nfr,nch,3))
r2=np.zeros((nmd,nT,nfr,nch,3))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{file}/{mdls[i]}_{Ts[j]}.npz')
        r1[i,j]=q['rc1']
        r2[i,j]=q['rc2']

"""Plot"""
# fig,ax=plt.subplots(nmd,nT,figsize=(20,12))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# labels=['Fix','Mid','Flex']
# # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# # nodes=[0.00,1/4,2/4,3/4,1.00]
# # cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# cmap=plt.cm.rainbow
# norm=colors.Normalize(vmin=0,vmax=nch)

# for i in range(nmd):
#     for j in range(nT):
#         ax[i,j].plot(r1[i,j,:,0,0],r2[i,j,:,0,0],'.',linewidth=wth)
#         ax[i,j].autoscale()
#         ax[i,j].minorticks_on()
#         ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
#         ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
#         ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
#         ax[i,j].spines['bottom'].set_linewidth(wth)
#         ax[i,j].spines['top'].set_linewidth(wth)
#         ax[i,j].spines['left'].set_linewidth(wth)
#         ax[i,j].spines['right'].set_linewidth(wth)
#         if i==nmd-1:
#             ax[i,j].set_xlabel('Raw X ($\\mathrm{\\AA}$) at %s K'%int(Ts[j]),fontsize=size)
#         if j==0:
#             ax[i,j].set_ylabel('Processed X ($\\mathrm{\\AA}$) of %s Model'%labels[i],fontsize=size)

# plt.tight_layout()
# plt.savefig(f'{figname}.png',format='png')
# plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

fig=plt.figure(figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=np.min(Ts),vmax=np.max(Ts))

print(np.mean(np.linalg.norm(r1[0,0,1:,0,:]-r1[0,0,:-1,0,:],axis=-1)**2))
dr=np.linalg.norm(r1[0,0,1:,0,:]-r1[0,0,:-1,0,:],axis=-1)**2
print(np.mean(dr[dr<10000]))
print(np.mean(np.linalg.norm(r2[0,0,1:,0,:]-r2[0,0,:-1,0,:],axis=-1)**2))
ax=plt.subplot(211)
# ax.plot(r2[0,0,:,0,0],r1[0,0,:,0,0],'.',linewidth=wth)
# ax.plot(r1[0,0,:,0,0],'-',linewidth=wth)
# ax.plot(r2[0,0,:,0,0]-(0.2*np.arange(nfr)),'-',linewidth=wth)
ax.plot(r2[0,0,:,0,0],'-',linewidth=wth,label='Proc')
ax.plot(r2[0,0,:,0,0]-r1[0,0,:,0,0],'-',linewidth=wth,label='Proc-Raw')
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
ax.set_xlabel('Frame',fontsize=size)
ax.set_ylabel('X ($\\mathrm{\\AA}$)',fontsize=size)
ax.legend(loc='best',fontsize=size)

ax=plt.subplot(212)
# ax.plot(r2[0,0,:,0,0],r1[0,0,:,0,0],'.',linewidth=wth)
# ax.plot(r1[0,0,:,0,0],'-',linewidth=wth)
# ax.plot(r2[0,0,:,0,0]-(0.2*np.arange(nfr)),'-',linewidth=wth)
ax.plot((r1[0,0,1:,0,0]-r1[0,0,:-1,0,0])[6450:6550],'-',linewidth=wth,alpha=0.5,label='Raw')
ax.plot((r2[0,0,1:,0,0]-r2[0,0,:-1,0,0])[6450:6550],'-',linewidth=wth,alpha=0.5,label='Proc')
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
ax.set_xlabel('Frame',fontsize=size)
ax.set_ylabel('dX ($\\mathrm{\\AA}$)',fontsize=size)
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()