"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
files=['../../all_atom/rg/rg.pyw',
       '../../hps/rg/rg.pyw',
       'fus_1.55sgm_1.0sgm/rg_itr20.pyw',
       'fus_1.5sgm_0.82sgm/rg_itr20.pyw',
       'fus_1.37sgm_0.6sgm/rg_itr20.pyw']
nmd=len(files)-2
figname='rg_compare'
nhi=200

"""Read Data and Calculate"""
rg=np.zeros((nmd+2,nhi))
p_rg=np.zeros((nmd+2,nhi))
rg_mean=np.zeros(nmd+2)
rg_std=np.zeros(nmd+2)
for i in range(nmd+2):
    q=pyw(files[i],'Probability')
    rg[i]=q[0]
    p_rg[i]=q[1]
    rg_mean[i]=pyw(files[i],'Mean')[0,0]
    rg_std[i]=pyw(files[i],'Standard')[0,0]

"""Plot"""
plt.figure(figsize=(10,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
labels=['AA','HPS','Qch','Mid','Flx']
colors=['k',(0.5,0.5,0.5),'r','g','b']

ax=plt.subplot(1,2,1)
for i in range(nmd+2):
    ax.plot(rg[i],p_rg[i],color=colors[i],linewidth=wth,label=labels[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(11,29)
ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
ax.set_ylabel('Probability Density ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
ax.legend(loc='best',fontsize=size)

ax=plt.subplot(1,2,2)
ax.fill_between([0,nmd+1],[rg_mean[0]-rg_std[0],rg_mean[0]-rg_std[0]],[rg_mean[0]+rg_std[0],rg_mean[0]+rg_std[0]],color='k',linewidth=0,alpha=0.2)
for i in range(nmd+2):
    ax.errorbar(i,rg_mean[i],yerr=rg_std[i],
                    marker='.',markersize=0,markeredgewidth=wth,
                    ecolor=colors[i],elinewidth=wth,capsize=2*wth)
    ax.plot(i,rg_mean[i],'.',color=colors[i],markersize=5*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nmd+2),labels)
ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
