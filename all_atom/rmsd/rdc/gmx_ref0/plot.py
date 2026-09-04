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
dir='rmsd'
figname=dir
ninit=16
inits=np.arange(ninit)
nrep=4
reps=np.arange(nrep)
time=1e-4

"""Read Data and Calculate"""
def xvg(file):
    q=[[],[]]
    f=open(file)
    lsf=f.readlines()
    for i in range(len(lsf)):
        if lsf[i][0] not in ['#','@']:
            for j in range(i,len(lsf)):
                ls_f=lsf[j].strip('\n').split()
                q[0]+=[float(ls_f[0])]
                q[1]+=[float(ls_f[1])]
            break
    q=np.array(q)
    return q

t=[]
d=[]
for i in range(ninit):
    t.append([])
    d.append([])
    for j in range(nrep):
        q=xvg(f'{dir}/{inits[i]}_{reps[j]}.xvg')
        t[-1]+=[(np.arange(len(q[0])))*time]
        d[-1]+=[q[1]*10]

"""Plot"""
fig,ax=plt.subplots(4,4,figsize=(20,15))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20

for i in range(4):
    for j in range(4):
        for k in range(nrep):
            ax[i,j].plot(t[i*4+j][k],d[i*4+j][k],linewidth=wth)
        ax[i,j].autoscale()
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        ax[i,j].set_xscale('log')
        ax[i,j].set_ylim(-0.5,33.5)
        if i==3:
            ax[i,j].set_xlabel('Time ($\\mathrm{\\mu{}s}$)',fontsize=size)
        if j==0:
            ax[i,j].set_ylabel('RMSD ($\\mathrm{\\AA}$)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
