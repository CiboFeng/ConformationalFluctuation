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
file1='rg_itr20.pyw'
file2='../fus_wt_mom1_type_1.37/rg_itr30.pyw'
figname=file1[:-4]+file2.split('/')[-1][2:-4]
sgm=3.82

"""Read Data and Calculate"""
q=pyw(file1,'Probability')
rg1=q[0]
p_rg1=q[1]
q=np.array(pyw(file2,'Probability'))
rg2=q[0]*sgm
p_rg2=q[1]/sgm
q=pyw(file1,'Mean')
rg1_mean=q[0]
q=np.array(pyw(file2,'Mean'))
rg2_mean=q[0]*sgm

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

fig,ax=plt.subplots(1,1,figsize=(6,5))

ax.plot(rg1,p_rg1,'r',linewidth=wth,label='$R_0=3.82$')
ax.plot([rg1_mean,rg1_mean],[np.min(p_rg1),np.max(p_rg1)],'r--',linewidth=wth,label='Mean(3.82)')
ax.plot(rg2,p_rg2,'b',linewidth=wth,label='$R_0=1.0$')
ax.plot([rg2_mean,rg2_mean],[np.min(p_rg2),np.max(p_rg2)],'b--',linewidth=wth,label='Mean(1.0)')
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
ax.set_xlim(10,30)
ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
ax.set_ylabel('Probability Density ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
