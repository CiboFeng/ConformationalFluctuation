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
files=['../../all_atom/dist/dist.npz',
       'dist/1.55sgm_1.0sgm_itr20.npz',
       'dist/1.5sgm_0.82sgm_itr20.npz',
       'dist/1.37sgm_0.6sgm_itr20.npz']
figname='dist'
nmd=len(files)-1
nat=163
nhi=50

"""Read Data and Calculate"""
d2=np.zeros((nmd+1,nat,nat,nhi))
p_d2=np.zeros((nmd+1,nat,nat,nhi))
d2_mean=np.zeros((nmd+1,nat,nat))
d2_std=np.zeros((nmd+1,nat,nat))
d1=np.zeros((nmd+1,nat,nhi))
p_d1=np.zeros((nmd+1,nat,nhi))
d1_mean=np.zeros((nmd+1,nat))
d1_std=np.zeros((nmd+1,nat))
for i in range(nmd+1):
    q=np.load(files[i])
    d2[i]=q['d2']
    p_d2[i]=q['p_d2']
    d2_mean[i]=q['d2_mean']
    d2_std[i]=q['d2_std']
    d1[i]=q['d1']
    p_d1[i]=q['p_d1']
    d1_mean[i]=q['d1_mean']
    d1_std[i]=q['d1_std']

"""Plot"""
fig,axs=plt.subplots(2,2,figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.82
cbarwth=0.03
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0/4,1/4,2/4,3/4,4/4]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=0,vmax=nat)
label0=['AA','Qch','Mid','Flx']
npair=10
Dseq=np.arange(nat-1,0,-round(nat/npair))
IJ=[]
for i in range(npair):
    # I=round(nat/2)-round(nat/2/npair)*i
    # J=nat-I-1
    # I=random.randint(0,nat)
    # J=random.randint(0,I-2)
    I=random.randint(0,nat-1-Dseq[i])
    J=I+Dseq[i]
    IJ+=[(I,J)]
    print(abs(I-J))

for i in range(nmd+1):
    j,k=divmod(i,2)
    ax=axs[j,k]
    for l,m in IJ:
        normc=norm(abs(l-m))
        rgb=cmap(normc)
        ax.plot(d2[i,l,m],p_d2[i,l,m],color=rgb,linewidth=wth)
        ax.plot([d2_mean[i,l,m],d2_mean[i,l,m]],[np.min(p_d2[i,l,m]),np.max(p_d2[i,l,m])],'--',color=rgb,linewidth=wth)
    ax.plot([],[],'--',color='k',linewidth=wth,label='Mean')
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
    ax.set_xlim(-4,84)
    if i==nmd:
        ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Probability Density\nof %s Model ($\\mathrm{\\AA}^{-1}$)'%label0[i],fontsize=size)
    if i==0:
        ax.legend(loc='best',fontsize=size)

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

fig,axs=plt.subplots(2,2,figsize=(15,10))
wth=2
size=20
lenmaj=15
lenmin=8
figrig=0.82
cbarwth=0.03
xtick=0.1
ytick=20
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=nat)
# color0=[(1.0,0.8,0.6),(0.8,0.6,1.0),(0.8,1.0,0.6)]
color0=[(0.50,0.50,0.50),(0.82,0.66,0.82),(0.56,0.90,0.94),(0.63,0.76,0.58)]
label0=['AA','Qch','Mid','Flx']
color_=['r','b']
def reshp(X):
    lenX=np.shape(X)[0]
    x=[]
    for i in range(lenX):
        for j in range(lenX-i):
            x+=[X[j,i+j]]
    x=np.array(x)
    return x
Dseq=[]
for i in range(nat):
    for j in range(nat-i):
        Dseq+=[i]
dseq=np.arange(nat)

for i in range(2):
    ax=axs[0,i]
    for j in range(nmd+1):
        ax.scatter(reshp(d2_mean[j]),reshp(d2_std[j]),c=color0[j],s=40*wth,label=label0[j])
        ax.scatter(reshp(d2_mean[j]),reshp(d2_std[j]),c=Dseq,s=10*wth,cmap=cmap,norm=norm)
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
    ax.set_xlabel('Distance Mean ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Distance Deviation ($\\mathrm{\\AA}$)',fontsize=size)
    if i==1:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([4,10,30])
        ax.set_xticklabels(['$4\\times{}10^0$','$10^1$','$3\\times{}10^1$'])
    else:
        ax.legend(loc='best',fontsize=size)

    ax=axs[1,i]
    ax2=ax.twinx()
    for j in range(nmd+1):
        ax.plot(dseq,d1_mean[j],linewidth=3*wth,color=color_[0])
        ax.plot(dseq,d1_mean[j],linewidth=2*wth,color=color0[j],label=label0[j])
        ax.plot(dseq,d1_std[j],'--',linewidth=3*wth,color=color_[1])
        ax.plot(dseq,d1_std[j],linewidth=2*wth,color=color0[j])
        ax2.plot(dseq,d1_mean[j],linewidth=3*wth,color=color_[0])
        ax2.plot(dseq,d1_mean[j],linewidth=2*wth,color=color0[j])
        ax2.plot(dseq,d1_std[j],'--',linewidth=3*wth,color=color_[1])
        ax2.plot(dseq,d1_std[j],linewidth=2*wth,color=color0[j])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=color_[0],labelcolor=color_[0])
    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=color_[0],labelcolor=color_[0])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(0)
    ax.spines['left'].set_color(color_[0])
    ax.yaxis.label.set_color(color_[0]) 
    ax.set_xlabel('Sequence Distance',fontsize=size)
    if i==0:
        ax.set_ylabel('Distance Mean ($\\mathrm{\\AA}$)',fontsize=size)
    else:
        ax.set_xscale('log')
        ax.legend(loc='best',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=color_[1],labelcolor=color_[1])
    ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=color_[1],labelcolor=color_[1])
    # ax2.xaxis.set_major_locator(MultipleLocator(xtick))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax2.yaxis.set_major_locator(MultipleLocator(ytick))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.spines['left'].set_linewidth(0)
    ax2.spines['right'].set_linewidth(wth)
    ax2.spines['right'].set_color(color_[1])
    ax2.yaxis.label.set_color(color_[1])
    if i==1:
        ax2.set_ylabel('Distance Deviation ($\\mathrm{\\AA}$)',fontsize=size)
        ax2.set_xscale('log')

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

plt.savefig(f'{figname}_mean_std.png',format='png')
plt.savefig(f'{figname}_mean_std.pdf',format='pdf')
plt.show()