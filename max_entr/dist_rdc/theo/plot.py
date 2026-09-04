"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
file='dist_rdc_3.npz'
figname=file[:-4]

"""Read Data and Calculate"""
q=np.load(file)
D=q['D']
D_repr=q['D_repr']
nfr,nat=np.shape(D)[:2]

D_diff=np.abs(D_repr-D)
ord=int(np.floor(np.log10(np.nanmax(np.abs(D_diff)))))
D_diff/=10**ord

show=0
D_show=(D_repr[show]-D[show])
ord_show=int(np.floor(np.log10(np.nanmax(np.abs(D_show)))))
D_show/=10**ord_show

"""Plot"""
fig,ax=plt.subplots(1,1,figsize=(6,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

ax.plot(D.reshape(-1),D_diff.reshape(-1),'.',linewidth=wth)
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
# ax.set_xlim(-4,84)
ax.set_xlabel('Raw Distance ($\\mathrm{\\AA}$)',fontsize=size)
ax.set_ylabel('Distance Difference ($10^{%s}\\mathrm{\\AA}$)'%ord,fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,ax=plt.subplots(1,1,figsize=(6,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
print(np.nanmin(D_show),np.nanmax(D_show))
xmin=-np.nanmin(D_show)/np.nanmax(np.abs(D_show))
xmax=np.nanmax(D_show)/np.nanmax(np.abs(D_show))
color=[(0,xmin,0),(0,0,0),(xmax,0,0)]
nodes=[0,-np.nanmin(D_show)/(np.nanmax(D_show)-np.nanmin(D_show)),1]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize()

ax.imshow(D_show,cmap=cmap,norm=norm)
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
ax.set_xlabel('Residue Index',fontsize=size)
ax.set_ylabel('Residue Index',fontsize=size)

fig.subplots_adjust(right=0.75)
plt.tight_layout(rect=[0,0,0.75,1])
top=ax.get_position().y1
bot=ax.get_position().y0
cbar_ax=fig.add_axes([0.78,bot,0.03,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Reproduce Error ($10^{%s}\\mathrm{\\AA}$)'%ord_show,fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_mat.png',format='png')
plt.savefig(f'{figname}_mat.pdf',format='pdf')
plt.show()