"""Import Modules"""
import random
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
from pyw import pyw

""""""
file1='top_show.npy'
dir2='../../chn_netw/netw'
figname='chn_netw'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nch=100
Pcut=50.0
size_lim=[1,300]
nhi=20
nq=11

""""""
P=np.load(file1)
for i in range(nmd):
    for j in range(nT):
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

""""""
x=np.zeros((nmd,nT,nhi,nq))
p=np.zeros((nmd,nT,nhi,nq))
x_mean=np.zeros((nmd,nT,nq))
x_std=np.zeros((nmd,nT,nq))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir2}/{mdls[i]}_{Ts[j]}.npz')
        x[i,j]=q['x']
        p[i,j]=q['p']
        x_mean[i,j]=q['x_mean']
        x_std[i,j]=q['x_std']
idx_sort=[0,2,1,3,4,5,6,7,8,9,10]
x=x[:,:,:,idx_sort]
p=p[:,:,:,idx_sort]
x_mean=x_mean[:,:,idx_sort]
x_std=x_std[:,:,idx_sort]

""""""
fig=plt.figure(figsize=(12,12))
gs=fig.add_gridspec(nrows=1,ncols=2,width_ratios=[1,1],wspace=0.4)
Gs=[gs[0].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.4),
    gs[1].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.4)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].nrows)] for i in range(len(Gs))]
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=-0.1
cbarwth=0.0
cbarspcwth=0.0
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
lenbar=8
xtick=[[0,4],[3,0.1]]
ytick=[[0,0.2],[1,10]]
ctick=1
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=np.min(bc[0,0]),vmax=np.max(bc[0,0]))
clrs=['r','g','b']
ann=[['(A)','(B)'],['(C)','(D)']]

ax=axs[0][0]
ax.annotate(ann[0][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
pos=nx.spring_layout(G[0,0],seed=42)
nodes=nx.draw_networkx_nodes(G[0,0],pos,ax=ax,node_size=node_size[i,j],node_color=bc[i,j],
                                cmap=cmap,vmin=np.min(bc[0,0]),vmax=np.max(bc[0,0]),alpha=0.9)
edges=nx.draw_networkx_edges(G[0,0],pos,ax=ax,width=wth,edge_color='gray',alpha=0.6)
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.axis('off')

divider=make_axes_locatable(ax)
cax=divider.append_axes('bottom',size='5%',pad=0.1)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cax,orientation='horizontal')
cbar.set_label('$g^\\mathrm{bc}_I$',fontdict=font)
cbar.ax.xaxis.set_major_locator(MultipleLocator(0.01))
cbar.ax.tick_params(direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar.outline.set_linewidth(wth)

col=[0,1,1]
row=[1,0,1]
idx=[1,5,10]
q=['\\langle{}N_I^\\mathrm{deg}\\rangle_I','\\mathrm{skw.}(g^\\mathrm{bc}_I)_I','c_\\mathrm{swn}']
for i in range(3):
    ax=axs[col[i]][row[i]]
    ax.annotate(ann[col[i]][row[i]],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.plot(x[j,0,:,idx[i]],p[j,0,:,idx[i]],color=clrs[j],linewidth=wth,label=Mdls[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(xtick[col[i]][row[i]]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[col[i]][row[i]]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$%s$'%q[i],fontdict=font)
    ax.set_ylabel('$f(%s)$'%q[i],fontdict=font)
    if i==1:
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

fig.canvas.draw()
tight_bbox=fig.get_tightbbox(fig.canvas.get_renderer())
alllft=tight_bbox.x0
allrig=tight_bbox.x1
allbot=tight_bbox.y0
alltop=tight_bbox.y1
figwth,fighei=fig.get_size_inches()
boxlft=np.inf
boxrig=-np.inf
boxbot=np.inf
boxtop=-np.inf
for ax in fig.axes:    
    pos=ax.get_position()
    left=pos.x0*figwth
    right=left+pos.width*figwth
    bottom=pos.y0*fighei
    top=bottom+pos.height*fighei
    boxlft=min(boxlft,left)
    boxrig=max(boxrig,right)
    boxbot=min(boxbot,bottom)
    boxtop=max(boxtop,top)
outlft=boxlft-alllft
outrig=allrig-boxrig
outbot=boxbot-allbot
outtop=alltop-boxtop

plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')