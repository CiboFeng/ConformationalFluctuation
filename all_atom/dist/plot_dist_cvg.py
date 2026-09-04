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
file='dist_cvg.npz'
figname=file[:-4]
nat=163
nhi=50

"""Read Data"""
data=np.load(file)
x1=data['x1']
px1=data['px1']
x2=data['x2']
px2=data['px2']

dseq=np.arange(nat)

"""Plot"""
nrow=3
ncol=3
fig,ax=plt.subplots(nrow,ncol,figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
lenbar=8
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=min(dseq),vmax=max(dseq))

K=[[5,10,20],[30,40,50],[60,70,80]]
for i in range(nrow):
    for j in range(ncol):
        k=K[i][j]
        l=nat-k-1
        normc=norm(dseq[abs(k-l)])
        rgb=cmap(normc)
        ax[i,j].plot(x1[k,l],px1[k,l],color='k',linewidth=2*wth)
        ax[i,j].plot(x1[k,l],px1[k,l],color=rgb,linewidth=wth)
        ax[i,j].plot(x2[k,l],px2[k,l],color=(0.5,0.5,0.5),linewidth=2*wth)
        ax[i,j].plot(x2[k,l],px2[k,l],color=rgb,linewidth=wth)
        if (i,j)==(0,0):
            ax[i,j].plot([],[],color='k',linewidth=2*wth,label='Set 1')
            ax[i,j].plot([],[],color=(0.5,0.5,0.5),linewidth=2*wth,label='Set 2')
            ax[i,j].legend(loc='best',fontsize=size)
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
        if i==nrow-1:
            ax[i,j].set_xlabel(f'Distance',fontsize=size)
        if j==0:
            ax[i,j].set_ylabel(f'Probability Density',fontsize=size)
        # ax[i,j].set_xscale('log')
        # ax[i,j].set_yscale('log')

fig.subplots_adjust(right=0.85)
plt.tight_layout(rect=[0,0,0.85,1])
top=ax[0,ncol-1].get_position().y1
bot=ax[nrow-1,ncol-1].get_position().y0
cbar_ax=fig.add_axes([0.88,bot,0.03,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()