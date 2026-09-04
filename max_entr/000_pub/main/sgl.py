"""Import Modules"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
files1=['../../../all_atom/rg/rg.pyw',
       '../../rg/fus_1.55sgm_1.0sgm/rg_itr20.pyw',
       '../../rg/fus_1.5sgm_0.82sgm/rg_itr20.pyw',
       '../../rg/fus_1.37sgm_0.6sgm/rg_itr20.pyw']
nmd=len(files1)-1
files2=['../../../all_atom/dist/dist.npz',
       '../../dist/dist/1.55sgm_1.0sgm_itr20.npz',
       '../../dist/dist/1.5sgm_0.82sgm_itr20.npz',
       '../../dist/dist/1.37sgm_0.6sgm_itr20.npz']
file3='../../coeff_cont_prob//coeff_20/sim_cont_prob_fus.pyw'
dir4='../../diff/msd'
dir5='../../autocor/autocor_dist'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
labels=['AA','Qch','Mid','Flx']
figname='sgl'
nat=163
nhi_rg=200
nfr=10000
dt=10e-6*10000 ### ns

"""Radius of Gyration"""
rg=np.zeros((nmd+1,nhi_rg))
p_rg=np.zeros((nmd+1,nhi_rg))
rg_mean=np.zeros(nmd+1)
rg_std=np.zeros(nmd+1)
for i in range(nmd+1):
    q=pyw(files1[i],'Probability')
    rg[i]=q[0]
    p_rg[i]=q[1]
    rg_mean[i]=pyw(files1[i],'Mean')[0,0]
    rg_std[i]=pyw(files1[i],'Standard')[0,0]

"""Distance"""
d_mean=np.zeros((nmd+1,nat))
d_std=np.zeros((nmd+1,nat))
for i in range(nmd+1):
    q=np.load(files2[i])
    d_mean[i]=q['d1_mean']
    d_std[i]=q['d1_std']

"""Contact Probability"""
P=np.zeros((nmd,nat,nat))
for i in range(nmd):
    P_tmp=pyw(file3.replace('//',f'/fus_{mdls[i]}/'),'Probability')[0]
    for j in range(nat):
        for k in range(nat):
            if abs(j-k)>=4:
                P[i,j,k]=P_tmp[j*nat+k]
P[1:]-=P[0:1]

"""Mean-squared Displacement"""
t_d=np.zeros((nmd,nfr,nat))
d=np.zeros((nmd,nfr,nat))
for i in range(nmd):
    q=pyw(f'{dir4}/{mdls[i]}_dyn.pyw','Mean-squared')
    for j in range(nfr):
        for k in range(nat):
            t_d[i,j,k]=q[0,j*nat+k]
            d[i,j,k]=q[2,j*nat+k]
dm=np.mean(d,axis=-1)
# dmax=np.mean(d[:,200:5000],axis=1)
# dmid=np.sum(d[:,:10]*np.arange(nfr-1,nfr-11,-1)[np.newaxis,:,np.newaxis],axis=1)/np.sum(np.arange(nfr-1,nfr-11,-1))
dmid=d[:,5]

"""Relaxation Time"""
t_ac=np.arange(nfr)*dt
ac=np.zeros((nmd,nfr,nat,nat))
for i in range(nmd):
    ac[i]=np.load(f'{dir5}/{mdls[i]}_dyn.npy')

t_relxe=np.zeros((nmd,nat,nat))
for i in range(nmd):
    for j in range(nat):
        for k in range(j):
            for l in range(nfr-1):
                if ac[i,l,j,k]>1/np.e and ac[i,l+1,j,k]<=1/np.e:
                    t_relxe[i,j,k]=((ac[i,l+1,j,k]-1/np.e)*t_ac[l]+(1/np.e-ac[i,l,j,k])*t_ac[l+1])/(ac[i,l+1,j,k]-ac[i,l,j,k])
                    t_relxe[i,k,j]=t_relxe[i,j,k]
                    break

ac_mean=np.zeros((nmd,nfr,nat))
t_relxe_mean=np.zeros((nmd,nat))
for i in range(nat):
    for j in range(nat-i):
        ac_mean[:,:,i]+=ac[:,:,j,i+j]
        t_relxe_mean[:,i]+=t_relxe[:,j,i+j]
    ac_mean[:,:,i]/=nat-i
    t_relxe_mean[:,i]/=nat-i
t_relxe_mean[:,0]=np.nan

"""Plot"""
fig=plt.figure(figsize=(25,25))
gs=fig.add_gridspec(nrows=3,ncols=1,height_ratios=[1,1.5,1],hspace=0.3)
Gs=[gs[0].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3),
    gs[1].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3),
    gs[2].subgridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3)]
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
xtick=100
ytick=100
ctick=5
clrs=['r','g','b']
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap1=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm1=colors.LogNorm(vmin=1e-3,vmax=1)
color=[(0,0,1),(1,1,1),(1,0,0)]
nodes=[0/2,3/13,13/13]
cmap2=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm2=colors.Normalize(vmin=-3e-2,vmax=1e-1)
ann=[['(D)','(E)','(F)'],['(G)','(H)','(I)'],['(J)','(K)','(L)']]

ax=axs[0][0]
ax.annotate(ann[0][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd+1):
    if i==0:
        rgb='k'
    else:
        rgb=clrs[i-1]
    ax.plot(rg[i],p_rg[i],color=rgb,linewidth=wth,label=labels[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(3))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(10.5,22)
ax.set_xlabel('$R\\mathrm{_g~(\\AA)}$',fontdict=font)
ax.set_ylabel('$f(R\\mathrm{_g)~(\\AA^{-1})}$',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ax=axs[0][1]
ax.annotate(ann[0][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
ax.fill_between([0,nmd],[rg_mean[0]-rg_std[0],rg_mean[0]-rg_std[0]],[rg_mean[0]+rg_std[0],rg_mean[0]+rg_std[0]],color='k',linewidth=0,alpha=0.2)
for i in range(nmd+1):
    if i==0:
        rgb='k'
    else:
        rgb=clrs[i-1]
    ax.errorbar(i,rg_mean[i],yerr=rg_std[i],
                    marker='.',markersize=0,markeredgewidth=wth,
                    ecolor=rgb,elinewidth=wth,capsize=2*wth)
    ax.plot(i,rg_mean[i],'.',color=rgb,markersize=5*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nmd+1),labels,rotation=30,fontdict=font)
ax.set_ylabel('$R\\mathrm{_g~(\\AA)}$',fontdict=font)

ax=axs[0][2]
ax.annotate(ann[0][2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(1,nmd+1):
    ax.plot(np.arange(nat),d_mean[i]/d_mean[0],linewidth=wth,color=clrs[i-1])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$|i-j|$',fontdict=font)
ax.set_ylabel('$\\langle{}d_{ij}\\rangle/\\langle{}d_{ij}^\\mathrm{AA}\\rangle$',fontdict=font)
# ax.set_xscale('log')

ax=axs[1][0]
ax.annotate(ann[1][0],xy=(0.0,1.25),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
img=ax.imshow(P[0],cmap=cmap1,norm=norm1)
ax.invert_yaxis()
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$i$',fontdict=font)
ax.set_ylabel('$j$',fontdict=font)
cbar=plt.colorbar(img,orientation='horizontal',location='top')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelfontfamily=font['family'],labelsize=font['size'])
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
cbar.outline.set_linewidth(wth)
cbar.set_label('$\\langle{}p_{ij}^\\mathrm{Qch}\\rangle$',fontdict=font)

for i in range(1,nmd):
    ax=axs[1][i]
    ax.annotate(ann[1][i],xy=(0.0,1.25),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    img=ax.imshow(P[i],cmap=cmap2,norm=norm2)
    ax.invert_yaxis()
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$i$',fontdict=font)
    ax.set_ylabel('$j$',fontdict=font)
    cbar=plt.colorbar(img,orientation='horizontal',location='top')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.set_ticks([-0.03,0,0.1])
    cbar.set_ticklabels(['-0.03','0','0.1'])
    cbar.outline.set_linewidth(wth)
    cbar.set_label('$\\langle{}p_{ij}^\\mathrm{%s}\\rangle-\\langle{}p_{ij}^\\mathrm{Qch}\\rangle$'%Mdls[i],fontdict=font)

ax=axs[2][0]
ax.annotate(ann[2][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(t_d[i,:,0],dm[i],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_minor_locator(AutoMinorLocator(0))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$\\Delta{}t$ (ns)',fontdict=font)
ax.set_ylabel('$\\langle\\Delta{}r^2\\rangle~(\\mathrm{\\AA^2)}$',fontdict=font)
ax.set_xscale('log')
ax.set_yscale('log')

ax=axs[2][1]
ax.annotate(ann[2][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(np.arange(nat),dmid[i],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$i$',fontdict=font)
ax.set_ylabel('$\\langle\\Delta{}r^2(t=0.5\\mathrm{ns})\\rangle~(\\mathrm{\\AA^2)}$',fontdict=font)

ax=axs[2][2]
ax.annotate(ann[2][2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(np.arange(nat),t_relxe_mean[i],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$|i-j|$',fontdict=font)
ax.set_ylabel('$\\tau(d_{ij})$ (ns)',fontdict=font)

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

fig=plt.figure(figsize=(25,7))
gs=fig.add_gridspec(nrows=1,ncols=3,width_ratios=[1,1,1],wspace=0.3)
axs=[fig.add_subplot(gs[i]) for i in range(gs.ncols)]
marglft=0.5
margrig=0.5
margbot=0.5
margtop=0.5
spcwth=0.2
spchei=-0.1
cbarwth=0.8
cbarspcwth=0.5
wth=3
font={'family':'Arial','size':30}
lenmaj=15
lenmin=8
color=[(1,0,0),(1,0.7,0.7),(1,1,1),(0.7,0.7,1),(0,0,1)]
nodes=[0/4,1/4,2/4,3/4,4/4]
cmap1=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=10,vmax=20)
ann=['(A)','(B)','(C)']
for i in range(3):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.set_axis_off()
figwth,fighei=fig.get_size_inches()
plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+1.0+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)
top=axs[-1].get_position().y1
bot=axs[-1].get_position().y0
cbar_ax=fig.add_axes([1.0-(margrig+1.0+cbarwth)/figwth,bot,cbarwth/figwth,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap1,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$R\\mathrm{_g~(\\AA)}$',fontdict=font)
cbar.ax.yaxis.set_major_locator(MultipleLocator(5))
cbar.ax.tick_params(direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
cbar.outline.set_linewidth(wth)
plt.savefig(f'{figname}_blk.png',format='png')
plt.savefig(f'{figname}_blk.pdf',format='pdf')

