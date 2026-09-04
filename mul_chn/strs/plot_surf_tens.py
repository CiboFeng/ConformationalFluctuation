"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
file='surf_tens.npz'
figname='surf_tens'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=20
nset=100
unit=0.101325 ### mN/m

"""Read Data and Calculate"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

q=np.load(file)
F=q['F']
p_F=q['p_F']
F_mean=q['F_mean']
F_std=q['F_std']

"""Plot"""
fig,axs=plt.subplots(2,3,figsize=(15,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

for i in range(nT):
    ax=axs[divmod(i,3)]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(F[j,i]*unit,p_F[j,i]/unit,color=rgb,linewidth=wth)
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
    # ax.set_xlim(9,73)
    if i==nT-1:
        ax.set_xlabel('Surface Tension (mN/m)',fontsize=size)
    ax.set_ylabel('Probability Density\nat %i K (m/mN)'%Ts[i],fontsize=size)
    # if i==0:
    #     ax.legend(loc='best',fontsize=size)

ax=axs[1,1]  
spc=3
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.errorbar(np.arange(nT)*(nmd+1+spc)+i,F_mean[i,:]*unit,yerr=F_std[i,:]*unit/np.sqrt(nset),
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nT)*(nmd+1+spc)+(nmd-1)/2,np.array(Ts).astype(int))
ax.set_xlabel('Temperature (K)',fontsize=size)
ax.set_ylabel('Surface Tension (mN/m)',fontsize=size)

axs[1,2].set_visible(False)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Structure Fluctuation',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_1.png',format='png')
plt.savefig(f'{figname}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(2,2,figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

for i in range(nmd):
    ax=axs[divmod(i,2)]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(F[i,j]*unit,p_F[i,j]/unit,color=rgb,linewidth=wth)
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
    # ax.set_xlim(9,73)
    if i==nmd-1:
        ax.set_xlabel('Surface Tension (mN/m)',fontsize=size)
    ax.set_ylabel('Probability Density\nof Flu.=%.2f (m/mN)'%rg_flu[i],fontsize=size)
    # if i==0:
    #     ax.legend(loc='best',fontsize=size)

ax=axs[1,1]  
spc=3
for i in range(nT):
    rgb=cmap(i/(nT-1))
    ax.errorbar(np.arange(nmd)*(nT+1+spc)+i,F_mean[:,i]*unit,yerr=F_std[:,i]*unit/np.sqrt(nset),
                     linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(9,73)
ax.set_xticks(np.arange(nmd)*(nT+1+spc)+(nT-1)/2,np.around(10**rg_flu,2))
ax.set_xlabel('Structure Fluctuation',fontsize=size)
ax.set_ylabel('Surface Tension (mN/m)',fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nT)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Temperature (K)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nT),1-1/(2*nT),nT))
cbar.set_ticklabels([int(Ts[i]) for i in range(nT)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()