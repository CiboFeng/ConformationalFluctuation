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
file='cont_prob_distr.pyw'
figname=file[:-4]
nat=74
nhist=50

"""Read Data and Calculate"""
x0=np.zeros((nat,nhist))
px0=np.zeros((nat,nhist))
y0=np.zeros((nat,nhist))
py0=np.zeros((nat,nhist))
x=np.zeros((nat,nhist))
px=np.zeros((nat,nhist))
y=np.zeros((nat,nhist))
py=np.zeros((nat,nhist))

for i in range(4,nat):
    q=pyw(file,f'{i}-th):')
    x0[i]=q[0]
    px0[i]=q[1]
    y0[i]=q[2]
    py0[i]=q[3]
    x[i]=q[4]
    px[i]=q[5]
    y[i]=q[6]
    py[i]=q[7]

d1=np.arange(4,nat)
xmin=min(np.min(y0),np.min(y))

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=min(d1),vmax=max(d1))

fig,ax=plt.subplots(2,2,figsize=(10,7))
for i in range(nat-4):
    normc=norm(d1[i])
    rgb=cmap(normc)
    ax[0,0].plot(x0[i],px0[i],color=rgb,linewidth=wth)
ax[0,0].autoscale()
ax[0,0].minorticks_on()
ax[0,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[0,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[0,0].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,0].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,0].spines['bottom'].set_linewidth(wth)
ax[0,0].spines['top'].set_linewidth(wth)
ax[0,0].spines['left'].set_linewidth(wth)
ax[0,0].spines['right'].set_linewidth(wth)
ax[0,0].set_xlim(-0.06,1.06)
ax[0,0].set_xlabel('Reference Contact Probability',fontsize=size)
ax[0,0].set_ylabel('Probability Density',fontsize=size)

for i in range(nat-4):
    normc=norm(d1[i])
    rgb=cmap(normc)
    ax[0,1].plot(y0[i],py0[i],color=rgb,linewidth=wth)
ax[0,1].autoscale()
ax[0,1].minorticks_on()
ax[0,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[0,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[0,1].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,1].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,1].spines['bottom'].set_linewidth(wth)
ax[0,1].spines['top'].set_linewidth(wth)
ax[0,1].spines['left'].set_linewidth(wth)
ax[0,1].spines['right'].set_linewidth(wth)
ax[0,1].set_xlim(xmin-1,1)
ax[0,1].set_ylim(-0.02,0.32)
ax[0,1].set_xlabel('Reference Logarithmic\nContact Probability',fontsize=size)
ax[0,1].set_ylabel('Probability Density',fontsize=size)

for i in range(nat-4):
    normc=norm(d1[i])
    rgb=cmap(normc)
    ax[1,0].plot(x[i],px[i],color=rgb,linewidth=wth)
ax[1,0].autoscale()
ax[1,0].minorticks_on()
ax[1,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[1,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[1,0].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,0].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,0].spines['bottom'].set_linewidth(wth)
ax[1,0].spines['top'].set_linewidth(wth)
ax[1,0].spines['left'].set_linewidth(wth)
ax[1,0].spines['right'].set_linewidth(wth)
ax[1,0].set_xlim(-0.06,1.06)
ax[1,0].set_xlabel('Fitted Contact Probability',fontsize=size)
ax[1,0].set_ylabel('Probability Density',fontsize=size)

for i in range(nat-4):
    normc=norm(d1[i])
    rgb=cmap(normc)
    ax[1,1].plot(y[i],py[i],color=rgb,linewidth=wth)
ax[1,1].autoscale()
ax[1,1].minorticks_on()
ax[1,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[1,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[1,1].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,1].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,1].spines['bottom'].set_linewidth(wth)
ax[1,1].spines['top'].set_linewidth(wth)
ax[1,1].spines['left'].set_linewidth(wth)
ax[1,1].spines['right'].set_linewidth(wth)
ax[1,1].set_xlim(xmin-1,1)
ax[1,1].set_ylim(-0.02,0.32)
ax[1,1].set_xlabel('Fitted Logarithmic\nContact Probability',fontsize=size)
ax[1,1].set_ylabel('Probability Density',fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)
plt.tight_layout(rect=[0,0,0.85,1])

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# # ytick=0.002
# color=[(1,0.8,0),(1,0,0),(0,0,0)]
# nodes=[0.00,4/9,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=min(cs),vmax=max(cs))
#
# fig,ax=plt.subplots(figsize=(10,5))
# for i in range(nc):
#     normc=norm(cs[i])
#     rgb=cmap(normc)
#     ax.plot(Q[gs.index(g),i,Ts.index(T)],F[gs.index(g),i,Ts.index(T)],color=rgb,linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# axis=plt.gca()
# axis.xaxis.set_minor_locator(AutoMinorLocator(2))
# axis.yaxis.set_minor_locator(AutoMinorLocator(2))
# axis.spines['bottom'].set_linewidth(wth)
# axis.spines['top'].set_linewidth(wth)
# axis.spines['left'].set_linewidth(wth)
# axis.spines['right'].set_linewidth(wth)
# ax.set_xlim(-0.06,1.06)
# ax.set_xlabel('$Q_{\\mathrm{CTP-Binding}}$',fontsize=size)
# ax.set_ylabel('Free Energy ($k_{\\mathrm{B}}T$)',fontsize=size)
# # ax.set_yscale('log')
# # ax.yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))
# sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
# sm.set_array([])
# cbar=plt.colorbar(sm)
# cbar.set_label('Concentration (mol/L)',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.outline.set_linewidth(wth)
# plt.tight_layout()
#
# plt.savefig(f'{figname}_{g}_{T}.png',format='png')
# plt.savefig(f'{figname}_{g}_{T}.pdf',format='pdf')
# plt.show()
#
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# # ytick=0.002
# color=[(1,0.8,0),(1,0,0),(0,0,0)]
# nodes=[0.00,4/9,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=min(gs),vmax=max(gs))
#
# fig,ax=plt.subplots(figsize=(10,5))
# for i in range(ng):
#     normc=norm(gs[i])
#     rgb=cmap(normc)
#     ax.plot(Q[i,cs.index(c),Ts.index(T)],F[i,cs.index(c),Ts.index(T)],color=rgb,linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# axis=plt.gca()
# axis.xaxis.set_minor_locator(AutoMinorLocator(2))
# axis.yaxis.set_minor_locator(AutoMinorLocator(2))
# axis.spines['bottom'].set_linewidth(wth)
# axis.spines['top'].set_linewidth(wth)
# axis.spines['left'].set_linewidth(wth)
# axis.spines['right'].set_linewidth(wth)
# ax.set_xlim(-0.06,1.06)
# ax.set_xlabel('$Q_{\\mathrm{CTP-Binding}}$',fontsize=size)
# ax.set_ylabel('Free Energy ($k_{\\mathrm{B}}T$)',fontsize=size)
# # ax.set_yscale('log')
# # ax.yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))
# sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
# sm.set_array([])
# cbar=plt.colorbar(sm)
# cbar.set_label('LJ Strength $\epsilon$ (kcal/mol)',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.outline.set_linewidth(wth)
# plt.tight_layout()
#
# plt.savefig(f'{figname}_{c}_{T}.png',format='png')
# plt.savefig(f'{figname}_{c}_{T}.pdf',format='pdf')
# plt.show()