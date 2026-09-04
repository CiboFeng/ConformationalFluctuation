"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
files=['../../all_atom/agl/agl.pyw',
       'agl/1.55sgm_1.0sgm_itr20.pyw',
       'agl/1.5sgm_0.82sgm_itr20.pyw',
       'agl/1.37sgm_0.6sgm_itr20.pyw']
nmd=len(files)-1
figname='agl'
nat=163
nhi=100

"""Read Data and Calculate"""
Tht=np.zeros((nmd+1,nat-2,nhi))
p_Tht=np.zeros((nmd+1,nat-2,nhi))
Tht_mean=np.zeros((nmd+1,nat-2))
Tht_std=np.zeros((nmd+1,nat-2))
Tht_all=np.zeros((nmd+1,nhi))
p_Tht_all=np.zeros((nmd+1,nhi))
Tht_mean_all=np.zeros(nmd+1)
Tht_std_all=np.zeros(nmd+1)
for i in range(nmd+1):
    q=pyw(files[i],'Each-residue Angle, Its Probability Density:')
    for j in range(nat-2):
        for k in range(nhi):
            Tht[i,j,k]=q[0,j*nhi+k]
            p_Tht[i,j,k]=q[1,j*nhi+k]
    Tht_mean[i]=pyw(files[i],'Mean of Each-residue Angle:')[0]
    Tht_std[i]=pyw(files[i],'Standard Deviation of Each-residue Angle:')[0]
    q=pyw(files[i],'All Angle, Its Probability Density:')
    Tht_all[i]=q[0]
    p_Tht_all[i]=q[1]
    Tht_mean_all[i]=pyw(files[i],'Mean of All Angle:')
    Tht_std_all[i]=pyw(files[i],'Standard Deviation of All Angle:')

"""Plot"""
fig,axs=plt.subplots(2,2,figsize=(10,8))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.82
cbarwth=0.03
xtick=0.1
ytick=20
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=np.max(p_Tht))
labels=['AA','Qch','Mid','Flx']

for i in range(nmd+1):
    j,k=divmod(i,2)
    ax=axs[j,k]
    ax2=ax.twiny()
    ax.imshow(p_Tht[i].T,cmap=cmap)
    ax.invert_yaxis()
    ax.plot(np.arange(nat-2),Tht_mean[i]/np.pi*nhi,color=(0,1,1),linewidth=wth)
    ax.plot([],[],label=labels[i])
    ax2.plot(p_Tht_all[i],Tht_all[i]/np.pi*nhi,color='b',linewidth=wth)
    ax2.plot([np.min(p_Tht_all[i]),np.max(p_Tht_all[i])],[Tht_mean_all[i]/np.pi*nhi]*2,color=(0.5,0.5,0.5),linestyle=':',linewidth=wth)
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
    ax.set_yticks([0,nhi/2,nhi],[0,90,180])
    if j==1:
        ax.set_xlabel('Residue Index',fontsize=size)
    if k==0:
        ax.set_ylabel('Angle ($\\circ$)',fontsize=size)
    ax.legend(loc='lower right',fontsize=size,handlelength=0.0,handletextpad=0.0)
    ax2.minorticks_on()
    ax2.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color='b',labelcolor='b')
    ax2.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color='b',labelcolor='b')
    # ax2.xaxis.set_major_locator(MultipleLocator(xtick))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax2.yaxis.set_major_locator(MultipleLocator(ytick))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['bottom'].set_linewidth(0)
    ax2.spines['top'].set_linewidth(wth)
    ax2.spines['top'].set_color('b')
    ax2.xaxis.label.set_color('b')
    if j==0:
        ax2.set_xlabel('Probability Density',fontsize=size)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Density',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
