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
file=[['fus_wt_mom1_type_1.55sgm/rg_itr0.pyw',
       'fus_wt_mom1_type_1.55sgm/rg_itr20.pyw'],
      ['fus_wt_mom1_type_1.37/rg_itr0.pyw',
       'fus_wt_mom1_type_1.37/rg_itr30.pyw']]
figname='rg_bg'
nhist=50
sgm=[[1.0,1.0],[3.82,3.82]]

"""Read Data and Calculate"""
rg=np.zeros((2,2,nhist))
p_rg=np.zeros((2,2,nhist))
rg_mean=np.zeros((2,2))
for i in range(2):
    for j in range(2):
        q=np.array(pyw(file[i][j],'Probability'))
        rg[i,j]=q[0]*sgm[i][j]
        p_rg[i,j]=q[1]/sgm[i][j]
        q=np.array(pyw(file[i][j],'Mean'))
        rg_mean[i,j]=q[0][0]*sgm[i][j]

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
line=[['r:','r-'],['b:','b-']]
label=[['Hetero-0','Hetero'],['Homo-0','Homo']]

fig,axs=plt.subplots(1,2,figsize=(12,5))
ax=axs[0]
for i in range(2):
    for j in range(2):
        ax.plot(rg[i,j],p_rg[i,j],line[i][j],linewidth=wth,label=label[i][j])
# ax.plot([rg0_mean,rg0_mean],[np.min(p_rg0),np.max(p_rg0)],'--',linewidth=wth,label='Mean')
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(11,30)
ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
ax.set_ylabel('Probability Density ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
ax.legend(loc='best',fontsize=size)

ax=axs[1]
ax.bar([0,1,2,3],rg_mean.reshape(-1),width=0.6,color='skyblue',edgecolor='black',linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
xticks=[0,1,2,3]
xlabels=np.array(label).reshape(-1)
ax.set_xticks(xticks,xlabels,rotation=45)
ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
