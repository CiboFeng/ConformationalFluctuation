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
coefffile1='coeff_19/coeff_fus_wt.pyw'
coefffile2='coeff_18_tdp43/coeff_tdp43.pyw'
figname='coeff_tdp43'
res1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
Nr=len(res1)
excl=[1,4,6,8,9,10,11,13,16,17,19]
Nr_=Nr-len(excl)

""""""
q=np.array(pyw(coefffile1,'Coefficient'))
nmo=len(q)
Coeff1=np.zeros((Nr,Nr,nmo))
for i in range(Nr):
    for j in range(Nr):
        for k in range(nmo):
            Coeff1[i,j,k]=q[k,i*Nr+j]

q=np.array(pyw(coefffile2,'Coefficient'))
nmo=len(q)
Coeff2=np.zeros((Nr,Nr,nmo))
for i in range(Nr):
    for j in range(Nr):
        for k in range(nmo):
            Coeff2[i,j,k]=q[k,i*Nr+j]

Coeff1_=np.zeros((Nr_,Nr_,nmo))
Coeff2_=np.zeros((Nr_,Nr_,nmo))
i_=0
for i in range(Nr):
    if i not in excl:
        j_=0
        for j in range(Nr):
            if j not in excl:
                Coeff1_[i_,j_,:]=Coeff1[i,j,:]
                Coeff2_[i_,j_,:]=Coeff2[i,j,:]
                j_+=1
        i_+=1

mat=np.zeros((2*nmo,Nr,Nr))
mat[:nmo]=Coeff1.transpose((2,0,1))
mat[nmo:]=Coeff2.transpose((2,0,1))
vec=np.zeros((2*nmo,Nr_**2))
vec[:nmo]=Coeff1_.transpose((2,0,1)).reshape(nmo,-1)
vec[nmo:]=Coeff2_.transpose((2,0,1)).reshape(nmo,-1)

CC=np.zeros((2*nmo,2*nmo))
for i in range(2*nmo):
    for j in range(2*nmo):
        CC[i,j]=np.corrcoef(vec[i],vec[j])[0,1]

""""""
fig,ax=plt.subplots(2*nmo,2*nmo+1,figsize=((2*nmo+1)*6,2*nmo*5))
wth=3
size=30
lenmaj=15
lenmin=8
lenbar=8
xtick=100
ytick=100

label=['$\\alpha_{%s}^\\mathrm{%s}$'%(i+1,'FUS-WT') for i in range(nmo)]+['$\\alpha_{%s}^\\mathrm{%s}$'%(i+1,'TDP43') for i in range(nmo)]

for i in range(2*nmo):
    color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
    nodes=[0,1/2,1]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.Normalize(vmin=np.min(mat[i]),vmax=np.max(mat[i]))
    img=ax[i,0].imshow(mat[i],cmap=cmap,norm=norm)
    ax[i,0].xaxis.set_minor_locator(MultipleLocator(xtick))
    ax[i,0].yaxis.set_minor_locator(MultipleLocator(ytick))
    ax[i,0].set_ylabel('Type',fontsize=size)
    if i==2*nmo-1:
        ax[i,0].set_xlabel('Type',fontsize=size)
    else:
        ax[i,0].set_xticks([],[])
    ax[i,0].invert_yaxis()
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
    for spine in ax[i,0].spines.values():
        spine.set_linewidth(wth)

for i in range(2*nmo):
    for j in range(2*nmo):
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
        if i==2*nmo-1:
            ax[i,j+1].set_xlabel(f'{label[j]}',fontsize=size)
        if j==0:
            ax[i,j+1].set_ylabel(f'{label[i]}',fontsize=size)
        ax[i,j+1].legend(loc='best',handlelength=0.0,handletextpad=0.0,fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()
