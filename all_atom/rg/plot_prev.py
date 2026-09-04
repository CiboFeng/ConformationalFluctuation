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
file_prev='prev/rg.pyw'
file='rg.pyw'
figname='rg_prev'
nat=163
nhist=50

"""Read Data and Calculate"""
q=pyw(file_prev,'Probability')
rg_prev=q[0]
p_rg_prev=q[1]
q=pyw(file_prev,'Mean')
rg_mean_prev=q[0]
q=pyw(file,'Probability')
rg=q[0]
p_rg=q[1]
q=pyw(file,'Mean')
rg_mean=q[0]

"""Plot"""
fig,ax=plt.subplots(1,1,figsize=(6,4))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20

ax.plot(rg_prev,p_rg_prev,'r',linewidth=wth)
ax.plot(rg,p_rg,'b',linewidth=wth)
ax.plot([rg_mean_prev,rg_mean_prev],[np.min(p_rg_prev),np.max(p_rg_prev)],'r--',linewidth=wth,label='Previous')
ax.plot([rg_mean,rg_mean],[np.min(p_rg),np.max(p_rg)],'b--',linewidth=wth,label='New')
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
ax.set_xlim(13,19)
ax.set_xlabel('Radius of Gyration',fontsize=size)
ax.set_ylabel('Probability Density',fontsize=size)
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
