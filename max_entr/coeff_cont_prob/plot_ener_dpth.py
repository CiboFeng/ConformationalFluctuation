"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import minimize_scalar
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
figname='ener_dpth'
rcre=[[1.37,0.6],[1.5,0.82],[1.55,1.0]]
nmd=len(files)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
ntp=len(res)
eps=0.2
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

"""Read Data"""
C=np.zeros((nmd,ntp,ntp))
for i in range(nmd):
    C_tmp=pyw(files[i],'Coefficient')[0]
    for j in range(ntp):
        for k in range(ntp):
            C[i,j,k]=C_tmp[j*ntp+k]

def lj(x,alf,rc,re,mu=1,eps=0.2):
    soft=alf/2*(1-np.tanh(mu*(x-rc)))
    ljlmd=4*eps*((re/x)**12-(re/x)**6)+eps
    ljlmdcut=ljlmd*(1-np.heaviside(x-2**(1/6)*re,0))
    return soft+ljlmdcut

em=np.zeros((nmd+1,nat,nat))
for i in range(nat):
    for j in range(i+1):
        em[0,i,j]=-eps*(lmd[seq[i]]+lmd[seq[j]])/2
        em[0,j,i]=em[0,i,j]
for i in range(nmd):
    for j in range(nat):
        for k in range(j+1):
            sgmij=(sgm[seq[j]]+sgm[seq[k]])/2
            u=minimize_scalar(lj,args=(C[i,seq[j],seq[k]],rcre[i][0]*sgmij,rcre[i][1]*sgmij)).fun
            if u>=0:
                u=C[i,seq[j],seq[k]]
            em[i+1,j,k]=u
            em[i+1,k,j]=u

"""Plot"""
fig,axs=plt.subplots(2,2,figsize=(12,9))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.82
cbarwth=0.03
xtick=0.1
ytick=20
xmin=-np.min(em)/np.max(np.abs(em))
xmax=np.max(em)/np.max(np.abs(em))
color=[(0,0,xmin),(1,1,1),(xmax,0,0)]
nodes=[0,-np.min(em)/(np.max(em)-np.min(em)),1]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=np.min(em),vmax=np.max(em))
labels=['HPS','Qch','Mid','Flx']
print(np.min(em),np.max(em))

for i in range(nmd+1):
    j,k=divmod(i,2)
    ax=axs[j,k]
    ax.imshow(em[i],cmap=cmap,norm=norm)
    ax.plot([],[],label=labels[i])
    ax.invert_yaxis()
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
    if j==1:
        ax.set_xlabel(f'Residue Index',fontsize=size)
    if k==0:
        ax.set_ylabel(f'Residue Index',fontsize=size)
    ax.legend(loc='upper left',fontsize=size,handlelength=0.0,handletextpad=0.0)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Energy Minimum (kcal/mol)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()