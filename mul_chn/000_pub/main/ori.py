"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw

"""Set Arguments"""
dir1='../../ori/ori'
dir2='../../dst/dst'
dir3='../../ori/ori_rel'
dir4='../../ori/mag'
dir5='../../res_distr/res_distr'
figname='ori'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi_z=250
nhi_tht=50
nhi_phi=50
nhi_mag=50
nat=163

"""Orientation"""
z=np.zeros((nmd,nT,nhi_z))
tht=np.zeros((nmd,nT,nhi_tht))
ptht=np.zeros((nmd,nT,nhi_z,nhi_tht))
tht_mean=np.zeros((nmd,nT,nhi_z))
tht_std=np.zeros((nmd,nT,nhi_z))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir1}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        tht[i,j]=q['tht']
        ptht[i,j]=q['p']
        tht_mean[i,j]=q['tht_mean']
        tht_std[i,j]=q['tht_std']
tht*=180/np.pi
ptht/=180/np.pi
tht_mean*=180/np.pi
tht_std*=180/np.pi

z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir2}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]

p_rand=np.pi/360*np.sin(np.pi*tht[0,0]/180)
ptht_slc=ptht/np.sum(ptht*180/nhi_tht,axis=-1,keepdims=True)
dptht=2*np.arccos(np.clip(np.sum(np.sqrt(ptht_slc)*np.sqrt(p_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_tht,axis=-1),-1.0,1.0))

"""Relative Orientation"""
z=np.zeros((nmd,nT,nhi_z))
phi=np.zeros((nmd,nT,nhi_phi))
pphi=np.zeros((nmd,nT,nhi_z,nhi_phi))
phi_mean=np.zeros((nmd,nT,nhi_z))
phi_std=np.zeros((nmd,nT,nhi_z))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir3}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        phi[i,j]=q['tht']
        pphi[i,j]=q['p']
        phi_mean[i,j]=q['tht_mean']
        phi_std[i,j]=q['tht_std']
phi*=180/np.pi
pphi/=180/np.pi
phi_mean*=180/np.pi
phi_std*=180/np.pi

# p_rand=1/180*np.ones(nhi_phi)
pphi_slc=pphi/np.sum(pphi*180/nhi_phi,axis=-1,keepdims=True)
dpphi=2*np.arccos(np.clip(np.sum(np.sqrt(pphi_slc)*np.sqrt(p_rand[np.newaxis,np.newaxis,np.newaxis])*180/nhi_phi,axis=-1),-1.0,1.0))

"""Magnetization"""
z=np.zeros((nmd,nT,nhi_z))
mag=np.zeros((nmd,nT,nhi_mag,3))
pmag=np.zeros((nmd,nT,nhi_z,nhi_mag,3))
mag_mean=np.zeros((nmd,nT,nhi_z,3))
mag_std=np.zeros((nmd,nT,nhi_z,3))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir4}/{mdls[i]}_{Ts[j]}.npz')
        z[i,j]=q['z']
        mag[i,j]=q['mag']
        pmag[i,j]=q['p']
        mag_mean[i,j]=q['mag_mean']
        mag_std[i,j]=q['mag_std']

"""Residue Distribution"""
z_seq=np.zeros((nmd,nT,nhi_z))
pz_seq=np.zeros((nmd,nT,nat,nhi_z))
zabs_seq_mean=np.zeros((nmd,nT,nat))
zabs_seq_std=np.zeros((nmd,nT,nat))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir5}/{mdls[i]}_{Ts[j]}.npz')
        z_seq[i,j]=q['Z']
        pz_seq[i,j]=q['pZ'].T
        zabs_seq_mean[i,j]=q['zabs_mean']
        zabs_seq_std[i,j]=q['zabs_std']

"""Plot"""
fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=1,ncols=2,width_ratios=[1,1],wspace=0.2)
Gs=[gs[0].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3),
    gs[1].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3)]
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
xtick=100
ytick=100
ctick=1
clrs=['r','g','b']
ann=[['(C)','(E)'],['(G)','(H)']]

ax=axs[0][0]
ax.annotate(ann[0][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(z[i,0]/((z2[i,0]-z1[i,0])/2+d[i,0]),dptht[i,0],color=clrs[i],linewidth=wth,label=Mdls[i])
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(-1,1)
ax.set_xlabel('$z_\\mathrm{norm}$',fontdict=font)
ax.set_ylabel('$d_\\mathrm{FR}(f(\\theta),f_\\mathrm{rand}(\\theta))$',fontdict=font)
# leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
# for txt in leg.get_texts():
#     txt.set_rotation(90) 
#     txt.set_rotation_mode('anchor')
#     txt.set_verticalalignment('top')
leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    txt.set_color(hnd.get_color())

ax=axs[0][1]
ax.annotate(ann[0][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(z[i,0]/((z2[i,0]-z1[i,0])/2+d[i,0]),dpphi[i,0],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.4))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(-1,1)
ax.set_yscale('log')
ax.set_xlabel('$z_\\mathrm{norm}$',fontdict=font)
ax.set_ylabel('$d_\\mathrm{FR}(f(\\phi),f_\\mathrm{rand}(\\phi))$',fontdict=font)

ax=axs[1][0]
ax.annotate(ann[1][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(z[i,0]/((z2[i,0]-z1[i,0])/2+d[i,0]),mag_mean[i,0,:,-1],color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlim(-1,1)
ax.set_xlabel('$z_\\mathrm{norm}$',fontdict=font)
ax.set_ylabel('$\\langle\\chi_z\\rangle$',fontdict=font)

ax=axs[1][1]
ax.annotate(ann[1][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.plot(zabs_seq_mean[i,0]/((z2[i,0]-z1[i,0])/2+d[i,0]),np.arange(nat),color=clrs[i],linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$\\langle|z_\\mathrm{norm}|\\rangle$',fontdict=font)
ax.set_ylabel('$i$',fontdict=font)

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

fig=plt.figure(figsize=(25,5))
gs=fig.add_gridspec(nrows=1,ncols=4,width_ratios=[1,1,1,1],wspace=0.4)
axs=[fig.add_subplot(gs[i]) for i in range(gs.ncols)]
font={'family':'Arial','size':30}
ann=['(A)','(B)','(D)','(F)']
for i in range(4):
    ax=axs[i]
    ax.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    if i==0:
        ax.text(0.0,0.5,'N terminal C terminal $\\theta \\phi z$',fontfamily=font['family'],fontsize=font['size'])
    ax.set_axis_off()
figwth,fighei=fig.get_size_inches()
plt.subplots_adjust(left=(marglft+outlft)/figwth,
                    right=1.0-(margrig+outrig+cbarwth+cbarspcwth)/figwth,
                    bottom=(margbot+outbot)/fighei,
                    top=1.0-(margtop+outtop)/fighei)
plt.savefig(f'{figname}_blk.png',format='png')
plt.savefig(f'{figname}_blk.pdf',format='pdf')