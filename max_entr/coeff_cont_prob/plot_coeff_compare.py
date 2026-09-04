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
files=['fus_1.55sgm_1.0sgm/coeff_20/coeff_fus.pyw',
       'fus_1.5sgm_0.82sgm/coeff_20/coeff_fus.pyw',
       'fus_1.37sgm_0.6sgm/coeff_20/coeff_fus.pyw']
figname='coeff_compare'
nmd=len(files)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
ntp=len(res)
excl=[1,4,6,8,9,10,11,13,17,19]
ntp_=ntp-len(excl)

"""Read Data"""
C=np.zeros((nmd,ntp,ntp))
for i in range(nmd):
    C_tmp=pyw(files[i],'Coefficient')[0]
    for j in range(ntp):
        for k in range(ntp):
            C[i,j,k]=C_tmp[j*ntp+k]
Lmd=(lmd+np.array(lmd).reshape(-1,1))/2

Lmd_=np.zeros((ntp_,ntp_))
C_=np.zeros((nmd,ntp_,ntp_))
i_=0
for i in range(ntp):
    if i not in excl:
        j_=0
        for j in range(ntp):
            if j not in excl:
                Lmd_[i_,j_]=Lmd[i,j]
                C_[:,i_,j_]=C[:,i,j]
                j_+=1
        i_+=1

CC=np.zeros(nmd)
for i in range(nmd):
    CC[i]=np.corrcoef(Lmd_.reshape(-1),C_[i].reshape(-1))[0,1]

"""Plot"""
plt.figure(figsize=(15,9))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
lenbar=8
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
label0=['Qch','Mid','Flx']

for i in range(nmd):
    ax=plt.subplot(2,3,i+1)
    ax.plot(C[divmod(i,nmd)[1]].reshape(-1),C[divmod(i+1,nmd)[1]].reshape(-1),'.',linewidth=wth)
    ax.plot([np.min(C[divmod(i,nmd)[1]]),np.max(C[divmod(i,nmd)[1]])],[np.min(C[divmod(i,nmd)[1]]),np.max(C[divmod(i,nmd)[1]])],linewidth=wth)
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
    ax.set_xlabel(f'Energy Coefficient (kcal/mol)\nof {label0[divmod(i,3)[1]]} Model',fontsize=size)
    ax.set_ylabel(f'Energy Coefficient (kcal/mol)\nof {label0[divmod(i+1,3)[1]]} Model',fontsize=size)

for i in range(nmd):
    ax=plt.subplot(2,3,i+4)
    ax.plot(Lmd_.reshape(-1),C_[i].reshape(-1),'.',linewidth=wth,label=f'C.C.={CC[i]:.2f}')
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
    ax.set_xlabel(f'$\\lambda$ in HPS Model',fontsize=size)
    ax.set_ylabel(f'Energy Coefficient (kcal/mol)\nof {label0[i]} Model',fontsize=size)
    ax.legend(loc='lower right',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()