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
file='rg.pyw'
figname=file[:-4]
nat=163
nhist=50

"""Read Data and Calculate"""
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

ax.plot(rg,p_rg,'r',linewidth=wth)
ax.plot([rg_mean,rg_mean],[np.min(p_rg),np.max(p_rg)],'--',linewidth=wth,label='Mean')
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
ax.set_xlabel('Radius of Gyration',fontsize=size)
ax.set_ylabel('Probability Density',fontsize=size)
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
