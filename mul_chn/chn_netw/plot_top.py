"""Import Modules"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

""""""
dir0='../../max_entr/rg/rg_heat'
dir='../cont/cont_slab'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
figname='top'
nch=100
fr=0
Pcut=50.0
size_lim=[1,100]

""""""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

P=np.zeros((nmd,nT,nch,nch))
for i in range(nmd):
    for j in range(nT):
        P[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')['Pchn'][fr]
        for k in range(nch):
            P[i,j,k,k]=0.0
P=np.heaviside(P-Pcut,0.0)
G=np.empty((nmd,nT),dtype=object)
deg=np.zeros((nmd,nT,nch))
bc=np.zeros((nmd,nT,nch))
for i in range(nmd):
    for j in range(nT):
        G[i,j]=nx.from_numpy_array(P[i,j])
        deg[i,j]=list(dict(G[i,j].degree()).values())
        betweenness=nx.betweenness_centrality(G[i,j],normalized=True)
        bc[i,j]=list(betweenness.values())
node_size=(size_lim[1]-size_lim[0])/(np.max(deg)-np.min(deg))*(deg-np.min(deg))+size_lim[0]

"""Plot"""
fig,axs=plt.subplots(nmd,nT,figsize=(20,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.87
cbarwth=0.03
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=np.min(bc),vmax=np.max(bc))

for i in range(nmd):
    for j in range(nT):
        ax=axs[i,j]
        pos=nx.spring_layout(G[i,j],seed=42)
        nodes=nx.draw_networkx_nodes(G[i,j],pos,ax=ax,node_size=node_size[i,j],node_color=bc[i,j],
                                     cmap=cmap,vmin=np.min(bc),vmax=np.max(bc),alpha=0.9)
        edges=nx.draw_networkx_edges(G[i,j],pos,ax=ax,width=1,edge_color='gray',alpha=0.6)
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        if i==nmd-1:
            ax.set_xlabel('%i K'%Ts[j],fontsize=size)
        if j==0:
            ax.set_ylabel('Flu.=%.2f'%10**rg_flu[i],fontsize=size)
        # ax.axis('off')

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Betweenness Centrality',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()