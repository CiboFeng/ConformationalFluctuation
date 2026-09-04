""""""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.insert(0,'D:\\Work\\Code\\functions')
from xtc import xtc_rd

""""""
file='fus_con.pro'
xtcfile=f'{file[:-4]}_sft.xtc'
datfile='fus.dat'
figname=file[:-4]
xtcstep=10000

""""""
timestep=[]
nnd=[]
nd=[]
f=open(file)
lsf=f.readlines()
for i in range(len(lsf)):
    if lsf[i]=='ITEM: TIMESTEP\n' and lsf[i+2]=='ITEM: NUMBER OF NODES\n':
        timestep+=[int(lsf[i+1].strip('\n'))]
    if lsf[i]=='ITEM: NUMBER OF NODES\n':
        nnd+=[int(lsf[i+1].strip('\n'))]
    if lsf[i]=='ITEM: NODES\n':
        nd.append([])
        for j in range(i+1,i+1+nnd[-1]):
            ls_f=lsf[j].strip('\n').split()
            nd[-1].append([float(ls_f[k]) for k in [2,3,4]])
nnd=np.array(nnd)
nd=np.array(nd)
nfr=np.shape(nd)[0]
ndz=[]
for i in range(nfr):
    ndz.append(sorted(list(set(nd[i,:,2]))))
ndz=np.array(ndz)

r,tp,b=xtc_rd(xtcfile,datfile)
r=r[:-1]
b=b[:-1]
nfr=np.shape(r)[0]
r-=b[:,np.newaxis,:,1]/2
zmax=np.max(r[:,:,2],axis=1)
zmin=np.min(r[:,:,2],axis=1)

""""""
fig=plt.figure(figsize=(8,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20

ax=plt.subplot(1,1,1)
for i in range(np.shape(ndz)[1]):
    ax.plot(timestep,ndz[:,i],'.',linewidth=wth)
ax.plot((np.arange(nfr)+1)*xtcstep,zmin,linewidth=wth,color='k')
ax.plot((np.arange(nfr)+1)*xtcstep,zmax,linewidth=wth,color='k',label='Matter')
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
ax.set_xscale('log')
ax.set_xlabel('Timestep',fontsize=size)
ax.set_ylabel('Boundary',fontsize=size)
ax.legend(loc='upper left',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
