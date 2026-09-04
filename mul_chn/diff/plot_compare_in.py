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
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
dir1='msd_raw_in'
dir2='msd_res_in'
figname='compare_in'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
dt=10e-6*10000 ### ns

"""Read Data"""
t=np.arange(nfr)*dt
d1=np.zeros((nmd,nT,nfr,nch,3))
d2=np.zeros((nmd,nT,nfr,nch,3))
for i in range(nmd):
    for j in range(nT):
        d1[i,j]=np.load(f'{dir1}/{mdls[i]}_{Ts[j]}.npz')['d2mean']
        d2[i,j]=np.load(f'{dir2}/{mdls[i]}_{Ts[j]}.npy')[:,:,81]

d1[d1==0.0]=np.nan
d1[:,:,0]=0.0
d2[d2==0.0]=np.nan
d2[:,:,0]=0.0
dm1=np.nanmean(d1,axis=-2)
dm1xy=np.sum(dm1[:,:,:,:2],axis=-1)
dm2=np.nanmean(d2,axis=-2)
dm2xy=np.sum(dm2[:,:,:,:2],axis=-1)

# def lastnan(X,last):
#     num_idx=np.where(~np.isnan(X))[0]
#     last_idx=num_idx[-last:] if len(num_idx)>=last else num_idx
#     X[last_idx]=np.nan
# for i in range(nmd):
#     for j in range(nT):
#         lastnan(dm1xy[i,j],100)
#         lastnan(dm2xy[i,j],100)

"""Plot"""
fig,ax=plt.subplots(nmd,nT,figsize=(20,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
labels=['Fix','Mid','Flex']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nch)

for i in range(nmd):
    for j in range(nT):
        ax[i,j].plot(t,dm1xy[i,j,:],linewidth=wth,label='Raw')
        ax[i,j].plot(t,dm2xy[i,j,:],linewidth=wth,label='Proc')
        ax[i,j].autoscale()
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        if i==nmd-1:
            ax[i,j].set_xlabel('Lag Time (ns) at %s K'%int(Ts[j]),fontsize=size)
        if j==0:
            ax[i,j].set_ylabel('MSD ($\\mathrm{\\AA}^2$) of %s Model'%labels[i],fontsize=size)
        ax[i,j].set_xscale('log')
        ax[i,j].set_yscale('log')
        # ax[i,j].set_ylim(40,1500)
        if (i,j)==(0,0):
            ax[i,j].legend(loc='upper left',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,ax=plt.subplots(nmd,nT,figsize=(20,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
labels=['Fix','Mid','Flex']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nch)

for i in range(nmd):
    for j in range(nT):
        ax[i,j].plot(t,dm1xy[i,j,:]/t/6,linewidth=wth,label='Raw')
        ax[i,j].plot(t,dm2xy[i,j,:]/t/6,linewidth=wth,label='Proc')
        ax[i,j].autoscale()
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        if i==nmd-1:
            ax[i,j].set_xlabel('Lag Time (ns) at %s K'%int(Ts[j]),fontsize=size)
        if j==0:
            ax[i,j].set_ylabel('Diffusion Coefficient\n($\\mathrm{\\AA}^2$/ns) of %s Model'%labels[i],fontsize=size)
        ax[i,j].set_xscale('log')
        # ax[i,j].set_yscale('log')
        if (i,j)==(0,0):
            ax[i,j].legend(loc='upper right',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()