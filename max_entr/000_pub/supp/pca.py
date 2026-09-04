"""Import Modules"""
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
pywfile='../../dist_rdc/pca.pyw'
npzfile='../../dist_rdc/pca_shp.npz'
labels=['AA','HPS','Qch','Mid','Flx']
nmd=len(labels)-2
figname='pca'

"""Read Data and Calculate"""
perct=pyw(pywfile,'Percentage').reshape(-1)
pc=[]*(nmd+2)
for i in range(nmd+2):
    pc+=[pyw(pywfile,f'({labels[i]})')]

data=np.load(npzfile)
cor=data['cor']
q_mean=data['q_mean']
q_std=data['q_std']
Q=data['Q']
P1=data['P1']
P2=data['P2']
nq=np.shape(cor)[1]

"""Plot"""
fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=3,ncols=1,height_ratios=[1,1,1],hspace=0.5)
Gs=[gs[0].subgridspec(nrows=1,ncols=5,width_ratios=[1,1,1,1,1],wspace=0.4),
    gs[1].subgridspec(nrows=1,ncols=4,width_ratios=[1,1,1,1],wspace=0.5),
    gs[2].subgridspec(nrows=1,ncols=5,width_ratios=[1,1,1,1,1],wspace=0.4)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=0.2
cbarwth=0.8
cbarspcwth=0.5
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
lenbar=8
xtick=[[],[1000,500],[]]
ytick=[[],[0.4,0.4,500,200],[]]
ctick=1
color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap1=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm1=colors.Normalize(vmin=0,vmax=np.max(P2[:,0,1]))
clrs=['k',(0.5,0.5,0.5),'r','g','b']
color=[(0,0,1),(1,1,1),(1,0,0)]
nodes=[0/2,1/2,2/2]
cmap2=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm2=colors.Normalize(vmin=-1,vmax=1)
ann=[['(A)','(B)','(C)','(D)','(E)'],['(F)','(G)','(H)','(I)'],['(J)','(K)','(L)','(M)','(N)']]
qs=['$\\zeta_1$','$\\zeta_2$','$R_1$','$R_2$','$R_3$','$R_\\mathrm{g}$','$\\mathit{\\Delta}$','$S$']

for i in range(nmd+2):
    ax=axs[0][i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    levels=np.linspace(0.01*np.max(P2[i,0,1])+0.99*np.min(P2[i,0,1]),np.max(P2[i,0,1]),20)
    ax.contourf(Q[i,0],Q[i,1],P2[i,0,1].T,levels=levels,cmap=cmap1)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(400))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(1.05*pc_lim[0,0]-0.05*pc_lim[0,1],1.05*pc_lim[0,1]-0.05*pc_lim[0,0])
    # ax.set_ylim(1.05*pc_lim[1,0]-0.05*pc_lim[1,1],1.05*pc_lim[1,1]-0.05*pc_lim[1,0])
    ax.set_xlim(-300,1300)
    ax.set_ylim(-500,500)
    ax.set_xlabel(f'$\\zeta_1$ (%i%%) of %s'%(perct[0]*100,labels[i]),fontdict=font)
    if i==0:
        ax.set_ylabel(f'$\\zeta_2$ (%i%%)'%(perct[1]*100),fontdict=font)

for i in range(2):
    ax=axs[1][i]
    ax.annotate(ann[1][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd+2):
        ax.plot(Q[j,i],P1[j,i]*100,color=clrs[j],linewidth=wth,label=labels[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(xtick[1][i]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\zeta_%s$'%(i+1),fontdict=font)
    ax.set_ylabel('$f(\\zeta_%s)\\times100$'%(i+1),fontdict=font)
    if i==0:
        ax.set_xlim(-300,1300)
        leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,ncol=2,columnspacing=1.0,prop={'family':font['family'],'size':font['size']*0.8})
        for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
            txt.set_color(hnd.get_color())

    ax=axs[1][i+2]
    ax.annotate(ann[1][i+2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.fill_between([0,nmd+1],[q_mean[0,i]-q_std[0,i],q_mean[0,i]-q_std[0,i]],[q_mean[0,i]+q_std[0,i],q_mean[0,i]+q_std[0,i]],color='k',linewidth=0,alpha=0.2)
    for j in range(nmd+2):
        ax.errorbar(j,q_mean[j,i],yerr=q_std[j,i],
                        marker='.',markersize=0,markeredgewidth=wth,
                        ecolor=clrs[j],elinewidth=wth,capsize=2*wth)
        ax.plot(j,q_mean[j,i],'.',color=clrs[j],markersize=5*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i+2]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(nmd+2),labels,rotation=45,fontdict=font)
    ax.set_ylabel('$\\zeta_%s$'%(i+1),fontdict=font)
    # ax.set_xlim(40,1500)

for i in range(nmd+2):
    ax=axs[2][i]
    ax.annotate(ann[2][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.imshow(cor[i],cmap=cmap2,norm=norm2)
    ax.invert_yaxis()
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size']*0.8)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size']*0.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    font_scl=font.copy()
    font_scl['size']=0.8*font['size']
    ax.set_xticks(np.arange(nq),qs,rotation=45,fontdict=font_scl)
    if i==0:
        ax.set_yticks(np.arange(nq),qs,fontdict=font_scl)
    else:
        ax.set_yticks([])
    ax.set_xlabel(labels[i],fontdict=font)

top=axs[0][-1].get_position().y1
bot=0.51
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap1)
sm.set_array([])
cbar1=plt.colorbar(sm,cax=cbar_ax)
cbar1.set_label('$f(\\zeta_1,\\zeta_2)$',fontdict=font)
cbar1.ax.set_yticks([])
cbar1.ax.tick_params(which='major',labelfontfamily=font['family'],labelsize=font['size'],direction='in',width=wth,length=lenmin)
cbar1.ax.tick_params(which='minor',labelfontfamily=font['family'],labelsize=font['size'],direction='in',width=wth,length=0)
cbar1.outline.set_linewidth(wth)

top=0.49
bot=axs[-1][-1].get_position().y0
cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap2,norm=norm2)
sm.set_array([])
cbar2=plt.colorbar(sm,cax=cbar_ax)
cbar2.set_label('CC.',fontdict=font)
cbar2.ax.yaxis.set_major_locator(MultipleLocator(1.0))
cbar2.ax.tick_params(labelfontfamily=font['family'],labelsize=font['size'],direction='in',width=wth,length=lenmin)
cbar2.outline.set_linewidth(wth)

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
pos1=cbar1.ax.get_position()
pos2=cbar2.ax.get_position()
boxrig=max(boxrig,(pos1.x0+pos1.width)*figwth,(pos2.x0+pos2.width)*figwth)
outlft=boxlft-alllft
outrig=allrig-boxrig
outbot=boxbot-allbot
outtop=alltop-boxtop

plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)
cbar1.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,
                       0.5-(margtop+outtop-margbot-outbot)/fighei/2+spchei/fighei/2,
                       cbarwth/figwth,
                       0.5-(margtop+outtop+margbot+outbot)/fighei/2-spchei/fighei/2])
cbar2.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,
                       (margbot+outbot)/fighei,
                       cbarwth/figwth,
                       0.5-(margtop+outtop+margbot+outbot)/fighei/2-spchei/fighei/2])

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

xytick=[[300,250,10,2,2,5,0.2,0.5],
        [1000,500,25,10,5,25,0.3,1],
        [250,200,5,2,2,5,0.2,0.5],
        [500,300,10,3,2,10,0.2,0.5],
        [500,400,10,5,5,10,0.3,1]]
for i in range(nmd+2):
    fig,axs=plt.subplots(nq,nq,figsize=(30,25))
    marglft=0.5
    margrig=0.5
    margbot=0.5
    margtop=0.5
    spcwth=0.1
    spchei=0.1
    cbarwth=0.8
    cbarspcwth=0.5
    wth=3
    font={'family':'Arial','size':30}
    lenmaj=15
    lenmin=8
    lenbar=8
    xtick=20
    ytick=20
    ctick=1
    color=[(1,1,1),(0,1,1),(0,0,1),(0,0,0)]
    nodes=[0/3,1/3,2/3,3/3]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.Normalize(vmin=0,vmax=np.max(P2))
    qs=['$\\zeta_1$','$\\zeta_2$','$R_1$','$R_2$','$R_3$','$R_\\mathrm{g}$','$\\mathit{\\Delta}$','$S$']
    us=['','',' $(\\mathrm{\\AA})$',' $(\\mathrm{\\AA})$',' $(\\mathrm{\\AA})$',' $(\\mathrm{\\AA})$','','']

    for j in range(nq):
        for k in range(nq):
            ax=axs[j,k]
            if j!=k:
                levels=np.linspace(0.01*np.max(P2[i,j,k])+0.99*np.min(P2[i,j,k]),np.max(P2[i,j,k]),20)
                ax.contourf(Q[i,k],Q[i,j],P2[i,k,j].T,levels=levels,cmap=cmap)
                ax.autoscale()
                ax.minorticks_on()
                ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
                ax.xaxis.set_major_locator(MultipleLocator(xytick[i][k]))
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))
                ax.yaxis.set_major_locator(MultipleLocator(xytick[i][j]))
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
                ax.spines['bottom'].set_linewidth(wth)
                ax.spines['top'].set_linewidth(wth)
                ax.spines['left'].set_linewidth(wth)
                ax.spines['right'].set_linewidth(wth)
                if j==nq-1:
                    ax.set_xlabel(qs[k]+us[k],fontdict=font)
                else:
                    ax.set_xticks([])
                if k==0:
                    ax.set_ylabel(qs[j]+us[j],fontdict=font)
                else:
                    ax.set_yticks([])
            else:
                ax.plot(Q[i,j],P1[i,j],linewidth=wth)
                ax.autoscale()
                ax.minorticks_on()
                ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
                ax.xaxis.set_major_locator(MultipleLocator(xytick[i][k]))
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))
                ax.yaxis.set_major_locator(MultipleLocator(xytick[i][j]))
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
                ax.spines['bottom'].set_linewidth(wth)
                ax.spines['top'].set_linewidth(wth)
                ax.spines['left'].set_linewidth(wth)
                ax.spines['right'].set_linewidth(wth)
                if j==nq-1:
                    ax.set_xlabel(qs[k]+us[k],fontdict=font)
                else:
                    ax.set_xticks([])
                if k==0:
                    ax.set_ylabel(qs[j]+us[j],fontdict=font)
                ax.set_yticks([])
    for j in range(nq):
        axs[j,j].set_xlim(axs[j-1,j].get_xlim())

    top=axs[0,-1].get_position().y1
    bot=axs[-1,-1].get_position().y0
    cbar_ax=fig.add_axes([1.0,bot,0.1,top-bot])
    sm=plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar=plt.colorbar(sm,cax=cbar_ax)
    cbar.set_label('$f(x,y)~|~x,y\\in${%s}'%(','.join(qs)),fontdict=font)
    cbar.ax.set_yticks([])
    cbar.ax.tick_params(direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.outline.set_linewidth(wth)

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
    pos=cbar.ax.get_position()
    boxrig=max(boxrig,(pos.x0+pos.width)*figwth)
    outlft=boxlft-alllft
    outrig=allrig-boxrig
    outbot=boxbot-allbot
    outtop=alltop-boxtop

    axs_2d=np.atleast_2d(axs)
    nrows,ncols=np.shape(axs_2d)
    plt.subplots_adjust(left=(marglft+outlft)/figwth,
                        right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                        bottom=(margbot+outbot)/fighei,
                        top=1.0-(margtop+outtop)/fighei,
                        hspace=nrows*spchei/(fighei-margbot-outbot-margtop-outtop-(nrows-1)*spchei),
                        wspace=ncols*spcwth/(figwth-marglft-outlft-margrig-outrig-cbarwth-cbarspcwth-(ncols-1)*spcwth))
    cbar.ax.set_position([1.0-(margrig+outrig+cbarwth)/figwth,(margbot+outbot)/fighei,cbarwth/figwth,1.0-(margtop+outtop+margbot+outbot)/fighei])

    plt.savefig(f'{figname}_{labels[i].lower()}.png',format='png',dpi=50)
    plt.savefig(f'{figname}_{labels[i].lower()}.pdf',format='pdf')
    # plt.show()