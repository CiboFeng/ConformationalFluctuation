"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
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
dir='../../coeff_cont_prob'
file='coeff_20/coeff_fus.pyw'
figname='coeff_cor'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
ntp=len(res)
res_=['A','N','D','Q','G','M','P','S','T','Y']
res_sort_=['D','S','T','N','Q','G','P','A','M','Y']
ntp_=len(res_)

"""Read Data"""
C=np.zeros((nmd,ntp,ntp))
for i in range(nmd):
    C_tmp=pyw(f'{dir}/fus_{mdls[i]}/{file}','Coefficient')[0]
    for j in range(ntp):
        for k in range(ntp):
            C[i,j,k]=C_tmp[j*ntp+k]
Lmd=(lmd+np.array(lmd).reshape(-1,1))/2

idx_sort=[res.index(i) for i in res_sort]
Lmd=Lmd[idx_sort]
Lmd=Lmd[:,idx_sort]
C=C[:,idx_sort]
C=C[:,:,idx_sort]

Lmd_=np.zeros((ntp_,ntp_))
C_=np.zeros((nmd,ntp_,ntp_))
i_=0
for i in range(ntp):
    if res_sort[i] in res_:
        j_=0
        for j in range(ntp):
            if res_sort[j] in res_:
                Lmd_[i_,j_]=Lmd[i,j]
                C_[:,i_,j_]=C[:,i,j]
                j_+=1
        i_+=1

CC=np.zeros(nmd)
for i in range(nmd):
    CC[i]=np.corrcoef(Lmd_.reshape(-1),C_[i].reshape(-1))[0,1]

"""Plot"""
# plt.figure(figsize=(10,2/3**0.5*10))
# marglft=0.5
# margrig=0.5
# margbot=0.5
# margtop=0.5
# spcwth=0.2
# spchei=-0.1
# cbarwth=0.0
# cbarspcwth=0.0
# wth=3
# font={'family':'Arial','size':30}
# lenmaj=15
# lenmin=8
# lenbar=8
# xtick=0.1
# ytick=20
# ctick=1

# ax=plt.subplot(1,1,1)
# vtx=[[0,1,1,0,-1,-1],[-1,-0.5,0.5,1,0.5,-0.5]]
# for i in range(6):
#     j=divmod(i,6)[1]
#     k=divmod(i+1,6)[1]
#     ax.plot([vtx[0][j],vtx[0][k]],[vtx[1][j],vtx[1][k]],color='k',linewidth=wth)
#     if divmod(i,2)[1]==0:
#         ax.plot([0,vtx[0][j]],[0,vtx[1][j]],color='k',linewidth=wth)

# # ax_lim=[]
# # for i in range(nmd):
# #     ax_lim.append([])
# #     ax_lim[i]+=[1.05*np.min(C[i])-0.05*np.max(C[i])]
# #     ax_lim[i]+=[-0.05*np.min(C[i])+1.05*np.max(C[i])]
# # ax_tick_maj=[[-0.5,0.0],[-0.5,0.0],[-0.5,0.0,0.5]]
# # ax_tick_min=[[-0.25],[-0.25,0.25],[-0.25,0.25]]
# ax_lim=[[1.05*np.min(C)-0.05*np.max(C),-0.05*np.min(C)+1.05*np.max(C)]]*3
# ax_tick_maj=[[-0.5,0.0,0.5]]*3
# ax_tick_min=[[-0.25,0.25]]*3

# def ax_tsf(x,y,n):
#     x=np.array(x)
#     y=np.array(y)
#     n1=divmod(n,3)[1]
#     n2=divmod(n+1,3)[1]
#     x_norm=(x-ax_lim[n1][0])/(ax_lim[n1][1]-ax_lim[n1][0])
#     y_norm=(y-ax_lim[n2][0])/(ax_lim[n2][1]-ax_lim[n2][0])
#     n3=divmod(2*n,6)[1]
#     n4=divmod(2*(n+1),6)[1]
#     X=x_norm*vtx[0][n3]+y_norm*vtx[0][n4]
#     Y=x_norm*vtx[1][n3]+y_norm*vtx[1][n4]
#     return X,Y

# for i in range(3):
#     for j in range(len(ax_tick_maj[i])):
#         xg10,yg10=ax_tsf(ax_tick_maj[i][j],ax_lim[divmod(i+1,3)[1]][0],i)
#         xg11,yg11=ax_tsf(ax_tick_maj[i][j],ax_lim[divmod(i+1,3)[1]][1],i)
#         xg20,yg20=ax_tsf(ax_lim[divmod(i-1,3)[1]][0],ax_tick_maj[i][j],i-1)
#         xg21,yg21=ax_tsf(ax_lim[divmod(i-1,3)[1]][1],ax_tick_maj[i][j],i-1)
#         ax.plot([xg10,xg11],[yg10,yg11],color='k',linewidth=wth,alpha=0.2)
#         ax.plot([xg20,xg21],[yg20,yg21],color='k',linewidth=wth,alpha=0.2)
# for i in range(3):
#     xo,yo=ax_tsf(np.mean(ax_lim[i]),1.05*ax_lim[divmod(i+1,3)[1]][0]-0.05*ax_lim[divmod(i+1,3)[1]][1],i)
#     ax.text(xo,yo,'$\\alpha_{IJ}$ of %s'%Mdls[i], 
#             va='top',ha='center',rotation=180/np.pi*np.arctan2(vtx[1][2*i],vtx[0][2*i]),rotation_mode='anchor',fontdict=font)
#     for j in range(len(ax_tick_maj[i])):
#         xt,yt=ax_tsf(ax_tick_maj[i][j],0.98*ax_lim[divmod(i+1,3)[1]][0]+0.02*ax_lim[divmod(i+1,3)[1]][1],i)
#         ax.text(xt,yt,ax_tick_maj[i][j], 
#                 va='bottom',ha='center',rotation=180/np.pi*np.arctan2(vtx[1][2*i],vtx[0][2*i]),rotation_mode='anchor',fontdict=font)

# for i in range(nmd):
#     x,y=ax_tsf(C[divmod(i,nmd)[1]].reshape(-1),C[divmod(i+1,nmd)[1]].reshape(-1),i)
#     ax.scatter(x,y,s=20*wth)
#     x0=[min(np.min(C[divmod(i,nmd)[1]]),np.min(C[divmod(i+1,nmd)[1]])),max(np.max(C[divmod(i,nmd)[1]]),np.max(C[divmod(i+1,nmd)[1]]))]
#     x,y=ax_tsf(x0,x0,i)
#     ax.plot(x,y,color=f'C{i}',linestyle='--',linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.spines['bottom'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.set_xticks([])
# ax.set_yticks([])
# # ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})


# plt.tight_layout()
# plt.savefig(f'{figname}.png',format='png',dpi=50)
# plt.savefig(f'{figname}.pdf',format='pdf')
# # plt.show()

fig=plt.figure(figsize=(25,18))
gs=fig.add_gridspec(nrows=2,ncols=1,height_ratios=[1,1.3],hspace=0.2)
Gs=[gs[0].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3),
    gs[1].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3)]
axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
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
xtick=[[],[0.05,0.2,0.05]]
ytick=0
ctick=1
color=[(0,0,1),(1,1,1),(1,0,0)]
nodes=[0/2,1/2,2/2]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=-0.6,vmax=0.6)
ann=[['(A)','(B)','(C)'],['(D)','(E)','(F)']]

for i in range(nmd):
    ax=axs[0][i]
    ax.annotate(ann[0][i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.scatter(C_[i].reshape(-1),C_[divmod(i+1,nmd)[1]].reshape(-1),s=20*wth)
    ax.plot([np.min(C_),np.max(C_)],[np.min(C_),np.max(C_)],'k--',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\alpha_{s_is_j}$ of %s'%Mdls[i],fontdict=font)
    ax.set_ylabel('$\\alpha_{s_is_j}$ of %s'%Mdls[divmod(i+1,nmd)[1]],fontdict=font)
    # leg=ax.legend(loc='upper right',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
    # for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    #     txt.set_color(hnd.get_color())

    ax=axs[1][i]
    ax.annotate(ann[1][i],xy=(0.0,1.2),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    img=ax.imshow(C_[divmod(i+1,nmd)[1]]-C_[i],cmap=cmap,norm=norm)
    ax.invert_yaxis()
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$s_i$',fontdict=font)
    ax.set_ylabel('$s_j$',fontdict=font)
    ax.set_xticks(np.arange(ntp_),res_sort_,fontdict=font)
    ax.set_yticks(np.arange(ntp_),res_sort_,fontdict=font)
    cbar=plt.colorbar(img,orientation='horizontal',location='top')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.outline.set_linewidth(wth)
    cbar.ax.xaxis.set_major_locator(MultipleLocator(0.5))
    cbar.set_label('$\\langle{}\\alpha_{s_is_j}^\\mathrm{%s}\\rangle-\\langle{}\\alpha_{s_is_j}^\\mathrm{%s}\\rangle$'%(Mdls[divmod(i+1,nmd)[1]],Mdls[i]),fontdict=font)

xlims=[axs[0][i].get_xlim() for i in range(nmd)]
ylims=[axs[0][i].get_ylim() for i in range(nmd)]
for i in range(nmd):
    axs[0][i].set_xlim(np.min(xlims),np.max(xlims))
    axs[0][i].set_ylim(np.min(ylims),np.max(ylims))

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

plt.savefig(f'{figname}.png',format='png',dpi=50)
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()