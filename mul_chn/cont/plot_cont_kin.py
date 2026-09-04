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

""""""
dir0='../../max_entr/rg/rg_heat'
file='cont_kin_show.npz'
figname='cont_kin'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)

""""""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

q=np.load(file)
Pintra1=q['Pintra1'][:,0]
Pinter1=q['Pinter1'][:,0]
Pintra2=q['Pintra2']
Pinter2=q['Pinter2']
Pintra2=np.concatenate((Pintra2,Pintra1[:,np.newaxis]),axis=1)
Pinter2=np.concatenate((Pinter2,Pinter1[:,np.newaxis]),axis=1)
nfr_=np.shape(Pintra2)[1]

"""Plot"""
fig,axs=plt.subplots(nmd,nfr_,figsize=(20,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=1)
frms=['$10^{-1}$ ns','$10^{0}$ ns','$10^{1}$ ns','$10^{2}$ ns','$10^{3}$ ns','Equ']

for i in range(nmd):
    for j in range(nfr_):
        ax=axs[i,j]
        ax.imshow(Pintra2[i,j],cmap=cmap,norm=norm)
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
        # ax.set_xscale('log')
        # ax.set_xlim(-4,84)
        if i==nmd-1:
            ax.set_xlabel(frms[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_intra.png',format='png')
plt.savefig(f'{figname}_intra.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nfr_,figsize=(20,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=5e-3)

for i in range(nmd):
    for j in range(nfr_):
        ax=axs[i,j]
        ax.imshow(Pinter2[i,j],cmap=cmap,norm=norm)
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
        # ax.set_xscale('log')
        # ax.set_xlim(-4,84)
        if i==nmd-1:
            ax.set_xlabel(frms[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Contact Probability',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(labelsize=size,which='minor',direction='in',width=wth,length=0)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}_inter.png',format='png')
plt.savefig(f'{figname}_inter.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(nmd,nfr_-1,figsize=(20,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-4,vmax=5e-3)

for i in range(nmd):
    for j in range(nfr_-1):
        ax=axs[i,j]
        ax.scatter(Pinter2[i,-1]*1e3,Pinter2[i,j]*1e3,linewidth=wth)
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
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        if i==nmd-1:
            ax.set_xlabel('Contact Probability\n$\\times{}10^{-3}$ (%s)'%frms[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Contact Probability\n$\\times{}10^{-3}$ of Flu.=%.2f'%10**rg_flu[i],fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_inter_cor.png',format='png')
plt.savefig(f'{figname}_inter_cor.pdf',format='pdf')
# plt.show()
