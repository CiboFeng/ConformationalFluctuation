import numpy as np
from scipy.spatial.distance import squareform,pdist
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
from matplotlib.pyplot import MultipleLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

""""""
coefffile='coeff_20/coeff_tdp43.pyw'
figname='coeff'
res1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
Nr=len(res1)
excl=[4,8,10,16,19]
Nr_=Nr-len(excl)

""""""
Sgm=np.zeros((Nr,Nr))
Lmd=np.zeros((Nr,Nr))
Chrg=np.zeros((Nr,Nr))
for i in range(Nr):
    for j in range(Nr):
        Sgm[i,j]=(sgm[i]+sgm[j])/2
        Lmd[i,j]=(lmd[i]+lmd[j])/2
        Chrg[i,j]=chrg[i]*chrg[j]
        
q=pyw(coefffile,'Coefficient')
nmo=len(q)
Coeff=np.zeros((Nr,Nr,nmo))
for i in range(Nr):
    for j in range(Nr):
        for k in range(nmo):
            Coeff[i,j,k]=q[k,i*Nr+j]

Sgm_=np.zeros((Nr_,Nr_))
Lmd_=np.zeros((Nr_,Nr_))
Chrg_=np.zeros((Nr_,Nr_))
Coeff_=np.zeros((Nr_,Nr_,nmo))
i_=0
for i in range(Nr):
    if i not in excl:
        j_=0
        for j in range(Nr):
            if j not in excl:
                Sgm_[i_,j_]=Sgm[i,j]
                Lmd_[i_,j_]=Lmd[i,j]
                Chrg_[i_,j_]=Chrg[i,j]
                Coeff_[i_,j_,:]=Coeff[i,j,:]
                j_+=1
        i_+=1

mat=np.zeros((3+nmo,Nr,Nr))
mat[0]=Sgm
mat[1]=Lmd
mat[2]=Chrg
mat[3:]=Coeff.transpose((2,0,1))
vec=np.zeros((3+nmo,Nr_**2))
vec[0]=Sgm_.reshape(-1)
vec[1]=Lmd_.reshape(-1)
vec[2]=Chrg_.reshape(-1)
vec[3:]=Coeff_.transpose((2,0,1)).reshape(nmo,-1)

CC=np.zeros((3+nmo,3+nmo))
for i in range(3+nmo):
    for j in range(3+nmo):
        CC[i,j]=np.corrcoef(vec[i],vec[j])[0,1]

""""""
fig,ax=plt.subplots(3+nmo,3+nmo+1,figsize=(25,18))
wth=3
size=30
lenmaj=15
lenmin=8
lenbar=8
xtick=100
ytick=100

label=['$\\sigma$','$\\lambda$','$q$']+['$\\alpha_%s$'%(i+1) for i in range(nmo)]

for i in range(3+nmo):
    color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
    nodes=[0,1/2,1]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.Normalize(vmin=np.min(mat[i]),vmax=np.max(mat[i]))
    img=ax[i,0].imshow(mat[i],cmap=cmap,norm=norm)
    ax[i,0].xaxis.set_minor_locator(MultipleLocator(xtick))
    ax[i,0].yaxis.set_minor_locator(MultipleLocator(ytick))
    ax[i,0].set_ylabel('Type',fontsize=size)
    if i==3+nmo-1:
        ax[i,0].set_xlabel('Type',fontsize=size)
    else:
        ax[i,0].set_xticks([],[])
    ax[i,0].invert_yaxis()
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
    for spine in ax[i,0].spines.values():
        spine.set_linewidth(wth)

for i in range(3+nmo):
    for j in range(3+nmo):
        ax[i,j+1].plot(vec[j],vec[i],'r.',linewidth=wth)
        ax[i,j+1].plot([],[],label=f'r={CC[i,j]:.2f}')
        ax[i,j+1].autoscale()
        ax[i,j+1].minorticks_on()
        ax[i,j+1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j+1].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax[i,j+1].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j+1].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j+1].spines['bottom'].set_linewidth(wth)
        ax[i,j+1].spines['top'].set_linewidth(wth)
        ax[i,j+1].spines['left'].set_linewidth(wth)
        ax[i,j+1].spines['right'].set_linewidth(wth)
        if i==3+nmo-1:
            ax[i,j+1].set_xlabel(f'{label[j]}',fontsize=size)
        if j==0:
            ax[i,j+1].set_ylabel(f'{label[i]}',fontsize=size)
        ax[i,j+1].legend(loc='best',handlelength=0.0,handletextpad=0.0,fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()
