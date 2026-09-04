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
file0='../../../all_atom/rg/rg.pyw'
file='rg_itr20.pyw'
figname=file[:-4]

"""Read Data and Calculate"""
q=pyw(file0,'Probability')
rg0=q[0]
p_rg0=q[1]
q=pyw(file,'Probability')
rg=q[0]
p_rg=q[1]
q=pyw(file0,'Mean')
rg0_mean=q[0]
q=pyw(file,'Mean')
rg_mean=q[0]

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

fig,ax=plt.subplots(2,1,figsize=(6,7))
ax[0].plot(rg0,p_rg0,'r',linewidth=wth)
ax[0].plot([rg0_mean,rg0_mean],[np.min(p_rg0),np.max(p_rg0)],'--',linewidth=wth,label='Mean')
ax[0].autoscale()
ax[0].minorticks_on()
ax[0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[0].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].spines['bottom'].set_linewidth(wth)
ax[0].spines['top'].set_linewidth(wth)
ax[0].spines['left'].set_linewidth(wth)
ax[0].spines['right'].set_linewidth(wth)
ax[0].set_xlim(10,30)
# ax[0].set_xlabel('Radius of Gyration',fontsize=size)
ax[0].set_ylabel('Probability Density\nof AA Model ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
ax[0].legend(loc='best',fontsize=size)

ax[1].plot(rg,p_rg,'r',linewidth=wth)
ax[1].plot([rg_mean,rg_mean],[np.min(p_rg),np.max(p_rg)],'--',linewidth=wth,label='Mean')
ax[1].autoscale()
ax[1].minorticks_on()
ax[1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[1].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1].spines['bottom'].set_linewidth(wth)
ax[1].spines['top'].set_linewidth(wth)
ax[1].spines['left'].set_linewidth(wth)
ax[1].spines['right'].set_linewidth(wth)
ax[1].set_xlim(10,30)
ax[1].set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
ax[1].set_ylabel('Probability Density\nof ME Model ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
ax[1].legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
