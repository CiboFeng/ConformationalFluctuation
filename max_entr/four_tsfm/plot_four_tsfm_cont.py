"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
dir='four_tsfm_cont'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nfr=10000
nat=163
nsms=19 ### Should be an odd number.

"""Read Data"""
A=np.zeros((nmd,round(nfr/2),nat,nat))
for i in range(nmd):
    q=np.load(f'{dir}/{mdls[i]}_dyn.npz')
    if i==0:
        f=q['freq']
    A[i]=q['A']

f=f[1:]
A=A[:,1:]
hlf=round((nsms-1)/2)
lenax=np.shape(A)[1]
A_=np.zeros_like(A,dtype=float)
for i in range(nsms):
    A_[:,hlf:-hlf]+=A[:,i:np.where(i-nsms+1!=0,i-nsms+1,lenax)]
A_[:,hlf:-hlf]/=nsms
for i in range(1,hlf):
    A__=np.zeros_like(A,dtype=float)
    for j in range(2*i+1):
        A__[:,i:-i]+=A[:,j:np.where(-2*i+j!=0,-2*i+j,lenax)]
    A__/=2*i+1
    A_[:,i]=A__[:,i]
    A_[:,-i-1]=A__[:,-i-1]
A_[:,0]=A[:,0]
A_[:,-1]=A[:,-1]
A=A_

"""Plot"""
fig,axs=plt.subplots(1,3,figsize=(15,4))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
npair=round(nat/10)
dseq=np.arange(nat-1,0,-round(nat/npair))
IJ=[]
for i in range(len(dseq)):
    I=random.randint(0,nat-1-dseq[i])
    J=I+dseq[i]
    IJ+=[(I,J)]

for i in range(nmd):
    ax=axs[i]
    for j,k in IJ:
        normc=norm(abs(j-k))
        rgb=cmap(normc)
        ax.plot(1/f,A[i,:,j,k],color=rgb,linewidth=wth)
    ax.plot([],[],label=labels[i])
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
    ax.set_xlabel('Period (ns)',fontsize=size)
    if i==0:
        ax.set_ylabel(f'Amplitude',fontsize=size)
    ax.set_xscale('log')
    ax.set_ylim(0,61)
    ax.legend(loc='upper left',fontsize=size,handlelength=0.0,handletextpad=0.0)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=ax[-1].get_position().y1
bot=ax[-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
