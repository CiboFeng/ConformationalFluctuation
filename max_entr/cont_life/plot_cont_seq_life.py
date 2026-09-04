"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir='cont_seq_life'
# dir='cont_seq_life_lin'
figname=dir
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nat=163
nhi=20

"""Read Data and Calculate"""
T0=np.zeros((nmd,nat,nhi))
p_T0=np.zeros((nmd,nat,nhi))
T0_mean=np.zeros((nmd,nat))
T0_std=np.zeros((nmd,nat))
T1=np.zeros((nmd,nat,nhi))
p_T1=np.zeros((nmd,nat,nhi))
T1_mean=np.zeros((nmd,nat))
T1_std=np.zeros((nmd,nat))
for i in range(nmd):
    q=np.load(f'{dir}/{mdls[i]}_dyn.npz')
    T0[i]=q['T0']
    p_T0[i]=q['p_T0']
    T0_mean[i]=q['T0_mean']
    T0_std[i]=q['T0_std']
    T1[i]=q['T1']
    p_T1[i]=q['p_T1']
    T1_mean[i]=q['T1_mean']
    T1_std[i]=q['T1_std']

"""Plot"""
fig,axs=plt.subplots(3,2,figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.88
cbarwth=0.03
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
labels=['Qch','Mid','Flx']

for i in range(nmd):
    ax=axs[i,0]
    for j in range(nat):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(T0[i,j],p_T0[i,j],color=rgb,linewidth=wth)
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
    ax.set_xticks([-1,0,1,2,3],['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
    if i==nmd-1:
        ax.set_xlabel('Life Time of 0 (ns)',fontsize=size)
    ax.set_ylabel('Probability Density\nof %s Model ($\\mathrm{\\AA}^{-1}$)'%labels[i],fontsize=size)

    ax=axs[i,1]
    for j in range(nat):
        normc=norm(j)
        rgb=cmap(normc)
        ax.plot(T1[i,j],p_T1[i,j],color=rgb,linewidth=wth)
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
    ax.set_xticks([-1,0,1,2,3],['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
    if i==nmd-1:
        ax.set_xlabel('Life Time of 1 (ns)',fontsize=size)
    # ax.set_ylabel('Probability Density\nof %s Model ($\\mathrm{\\AA}^{-1}$)'%label0[i],fontsize=size)

fig.subplots_adjust(right=figrig)
plt.tight_layout(rect=[0,0,figrig,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
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

fig,axs=plt.subplots(1,2,figsize=(15,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
labels=['Qch','Mid','Flx']
colors=[(0.82,0.66,0.82),(0.56,0.90,0.94),(0.63,0.76,0.58)]
colors_=['r','b']

ax=axs[0]
ax2=ax.twinx()
for i in range(nmd):
    ax.plot(np.arange(2,nat),T0_mean[i,2:],color=colors_[0],linewidth=3*wth)
    ax.plot(np.arange(2,nat),T0_mean[i,2:],color=colors[i],linewidth=2*wth,label=labels[i])
    ax2.plot(np.arange(2,nat),T0_std[i,2:],color=colors_[1],linewidth=3*wth)
    ax2.plot(np.arange(2,nat),T0_std[i,2:],color=colors[i],linewidth=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=colors_[0],labelcolor=colors_[0])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=colors_[0],labelcolor=colors_[0])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(0)
ax.spines['left'].set_color(colors_[0])
ax.yaxis.label.set_color(colors_[0]) 
ax.set_xscale('log')
ax.set_ylim(-1.01,0.74)
# ax.set_ylim(-3,77)
ax.set_yticks([-1.0,-0.5,0.0,0.5],['$10^{-1.0}$','$10^{-0.5}$','$10^{0.0}$','$10^{0.5}$'])
ax.set_xlabel('Sequence Distance',fontsize=size)
ax.set_ylabel('Life Time Mean of 0 (ns)',fontsize=size)
ax.legend(loc='lower right',fontsize=size)
ax2.minorticks_on()
ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=colors_[1],labelcolor=colors_[1])
ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=colors_[1],labelcolor=colors_[1])
# ax2.xaxis.set_major_locator(MultipleLocator(xtick))
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.yaxis.set_major_locator(MultipleLocator(ytick))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(wth)
ax2.spines['right'].set_color(colors_[1])
ax2.yaxis.label.set_color(colors_[1])
ax2.set_ylim(-1.01,0.74)
# ax.set_ylim(-3,77)
ax2.set_yticks([-1.0,-0.5,0.0,0.5],['$10^{-1.0}$','$10^{-0.5}$','$10^{0.0}$','$10^{0.5}$'])
ax2.set_ylabel('Life Time Deviation of 0 (ns)',fontsize=size)

ax=axs[1]
ax2=ax.twinx()
for i in range(nmd):
    ax.plot(np.arange(2,nat),T1_mean[i,2:],color=colors_[0],linewidth=3*wth)
    ax.plot(np.arange(2,nat),T1_mean[i,2:],color=colors[i],linewidth=2*wth,label=labels[i])
    ax.plot(np.arange(2,nat),T1_std[i,2:],color=colors_[1],linewidth=3*wth)
    ax.plot(np.arange(2,nat),T1_std[i,2:],color=colors[i],linewidth=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=colors_[0],labelcolor=colors_[0])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=colors_[0],labelcolor=colors_[0])
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(0)
ax.spines['left'].set_color(colors_[0])
ax.yaxis.label.set_color(colors_[0]) 
ax.set_xscale('log')
ax.set_ylim(-1.01,0.74)
ax.set_yticks([-1.0,-0.5,0.0,0.5],['$10^{-1.0}$','$10^{-0.5}$','$10^{0.0}$','$10^{0.5}$'])
ax.set_xlabel('Sequence Distance',fontsize=size)
ax.set_ylabel('Life Time Mean of 1 (ns)',fontsize=size)
ax2.minorticks_on()
ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=colors_[1],labelcolor=colors_[1])
ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=colors_[1],labelcolor=colors_[1])
# ax2.xaxis.set_major_locator(MultipleLocator(xtick))
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.yaxis.set_major_locator(MultipleLocator(ytick))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.spines['left'].set_linewidth(0)
ax2.spines['right'].set_linewidth(wth)
ax2.spines['right'].set_color(colors_[1])
ax2.yaxis.label.set_color(colors_[1])
ax2.set_ylim(-1.01,0.74)
ax2.set_yticks([-1.0,-0.5,0.0,0.5],['$10^{-1.0}$','$10^{-0.5}$','$10^{0.0}$','$10^{0.5}$'])
ax2.set_ylabel('Life Time Deviation of 1 (ns)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_mean_std.png',format='png')
plt.savefig(f'{figname}_mean_std.pdf',format='pdf')
plt.show()