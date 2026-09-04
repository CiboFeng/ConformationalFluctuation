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
file='mom.pyw'
figname=file[:-4]
nat=74
nhist=50
nmo=3

"""Read Data and Calculate"""
q=np.array(pyw(file,'Distance,'))
d0=q[1:1+nmo]
p0=q[1+nmo:1+2*nmo]
d=q[1+2*nmo:1+3*nmo]
p=q[1+3*nmo:1+4*nmo]

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

fig,ax=plt.subplots(2,nmo,figsize=(4*nmo+3,7))
for i in range(nmo):
    ax[0,i].plot(d0[i],linewidth=wth,label='Reference')
    ax[0,i].plot(d[i],linewidth=wth,label='Fitted')
    ax[0,i].autoscale()
    ax[0,i].minorticks_on()
    ax[0,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[0,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax[0,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0,i].spines['bottom'].set_linewidth(wth)
    ax[0,i].spines['top'].set_linewidth(wth)
    ax[0,i].spines['left'].set_linewidth(wth)
    ax[0,i].spines['right'].set_linewidth(wth)
    # ax[0,i].set_xlim(-7,165)
    ax[0,i].set_xscale('log')
    ax[0,i].set_yscale('log')
    # ax[0,i].set_xlabel('Sequence Distance',fontsize=size)
    ax[0,i].set_ylabel(f'The {i+1}-th Moment\nof Distance',fontsize=size)
    if i==0:
        ax[0,i].legend(loc='best',fontsize=size)

    ax[1,i].plot(p0[i],linewidth=wth,label='Reference')
    ax[1,i].plot(p[i],linewidth=wth,label='Fitted')
    ax[1,i].autoscale()
    ax[1,i].minorticks_on()
    ax[1,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[1,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax[1,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1,i].spines['bottom'].set_linewidth(wth)
    ax[1,i].spines['top'].set_linewidth(wth)
    ax[1,i].spines['left'].set_linewidth(wth)
    ax[1,i].spines['right'].set_linewidth(wth)
    # ax[1,i].set_xlim(0.4,2.3)
    ax[1,i].set_xscale('log')
    ax[1,i].set_yscale('log')
    ax[1,i].set_xlabel('Sequence Distance',fontsize=size)
    ax[1,i].set_ylabel(f'The {i+1}-th Moment of\nContace Probability',fontsize=size)
    # ax[1,i].legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
