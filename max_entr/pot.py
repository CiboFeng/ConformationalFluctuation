""""""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

""""""
figname='pot'
sgm=5.88

""""""
def lj(x,eps=0.2,alf=-1.0,re=0.6*sgm,rc=1.0*sgm,mu=1):
    ljlmd=4*eps*((re/x)**12-(re/x)**6)+eps
    ljlmdcut=ljlmd*(1-np.heaviside(x-2**(1/6)*re,0))
    soft=alf/2*(1-np.tanh(mu*(x-rc)))
    return ljlmdcut+soft

alf=np.arange(-2,2,0.5)
rc=np.arange(1.0,1.9,0.2)*sgm
re=np.arange(0.6,1.1,0.1)*sgm

""""""
plt.figure(figsize=(15,5))
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
norm=colors.Normalize(vmin=np.min(alf),vmax=np.max(alf))
ax=plt.subplot(1,3,1)
x=np.arange(0.1,10,0.01)
for i in range(len(alf)):
    normc=norm(alf[i])
    rgb=cmap(normc)
    ax.plot(x,lj(x,alf=alf[i]),color=rgb,linewidth=wth)
ax.plot([],[],label='$r_\\mathrm{e}=%.1f\\mathrm{\\AA}$'%(0.6*sgm))
ax.plot([],[],label='$r_\\mathrm{c}=%.1f\\mathrm{\\AA}$'%(1.0*sgm))
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
ax.set_ylim(-2.4,3.4)
ax.set_xlabel('$r~(\\mathrm{\\AA})$',fontsize=size)
ax.set_ylabel('$U(r)$ (kcal/mol)',fontsize=size)
ax.legend(loc='upper right',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,ax=plt.gca())
cbar.set_label('$\\alpha$ (kcal/mol)',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

color=[(1,1,0),(1,0,0),(0,0,0)]
nodes=[0.00,1/2,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=np.min(rc),vmax=np.max(rc))
ax=plt.subplot(1,3,2)
x=np.arange(0.1,15,0.01)
for i in range(len(rc)):
    normc=norm(rc[i])
    rgb=cmap(normc)
    ax.plot(x,lj(x,rc=rc[i]),color=rgb,linewidth=wth)
ax.plot([],[],label='$\\alpha=%.1f\\mathrm{kcal/mol}$'%(-1.0))
ax.plot([],[],label='$r_\\mathrm{e}=%.1f\\mathrm{\\AA}$'%(0.6*sgm))
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
ax.set_ylim(-2.4,3.4)
ax.set_xlabel('$r~(\\mathrm{\\AA})$',fontsize=size)
ax.legend(loc='upper right',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,ax=plt.gca())
cbar.set_label('$r_\\mathrm{c}~(\\mathrm{\\AA})$',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

color=[(0,1,1),(0,0,1),(0,0,0)]
nodes=[0.00,1/2,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=np.min(re),vmax=np.max(re))
ax=plt.subplot(1,3,3)
x=np.arange(0.1,10,0.01)
for i in range(len(re)):
    normc=norm(re[i])
    rgb=cmap(normc)
    ax.plot(x,lj(x,re=re[i]),color=rgb,linewidth=wth)
ax.plot([],[],label='$\\alpha=%.1f\\mathrm{kcal/mol}$'%(-1.0))
ax.plot([],[],label='$r_\\mathrm{c}=%.1f\\mathrm{\\AA}$'%(1.0*sgm))
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
ax.set_ylim(-2.4,3.4)
ax.set_xlabel('$r~(\\mathrm{\\AA})$',fontsize=size)
ax.legend(loc='upper right',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,ax=plt.gca())
cbar.set_label('$r_\\mathrm{e}~(\\mathrm{\\AA})$',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

plt.figure(figsize=(8,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=plt.subplot(1,1,1)
x=np.arange(0.1,5,0.01)
ax.plot(x,lj(x,rc=1.37,re=0.6),linewidth=wth,color=cmap(1.0),label='$r_\\mathrm{c}=1.37\\sigma_{ij}$,$r_\\mathrm{e}=0.6\\sigma_{ij}$')
ax.plot(x,lj(x,rc=1.5,re=0.82),linewidth=wth,color=cmap(0.5),label='$r_\\mathrm{c}=1.5\\sigma_{ij}$,$r_\\mathrm{e}=0.82\\sigma_{ij}$')
ax.plot(x,lj(x,rc=1.55,re=1.0),linewidth=wth,color=cmap(0.0),label='$r_\\mathrm{c}=1.55\\sigma_{ij}$,$r_\\mathrm{e}=1.0\\sigma_{ij}$')
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
ax.set_ylim(-2.4,3.4)
ax.set_xlabel('$r~(\\sigma_{ij})$',fontsize=size)
ax.set_ylabel('$U(r)$ (kcal/mol)',fontsize=size)
ax.legend(loc='upper right',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
plt.show()

def ljlmdcut(x,eps=0.2,re=0.6*sgm):
    ljlmd=4*eps*((re/x)**12-(re/x)**6)+eps
    ljlmdcut=ljlmd*(1-np.heaviside(x-2**(1/6)*re,0))
    return ljlmdcut
def soft(x,alf=-1.0,rc=1.0*sgm,mu=1):
    soft=alf/2*(1-np.tanh(mu*(x-rc)))
    return soft

plt.figure(figsize=(8,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20

ax=plt.subplot(1,1,1)
x=np.arange(0.1,5,0.01)
ax.plot(x,ljlmdcut(x,re=0.6),linewidth=wth,label='$U(r)=4\\epsilon[(\\frac{r_\\mathrm{e}}{r})^{12}-(\\frac{r_\\mathrm{e}}{r})^6]+\\epsilon$')
ax.plot(x,soft(x,rc=1.5),linewidth=wth,label='$U(r)=\\alpha/2[1-\\tanh\\mu(r-r_\\mathrm{c})]$')
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
ax.set_ylim(-1.4,3.4)
ax.set_xlabel('$r~(\\sigma_{ij})$',fontsize=size)
ax.set_ylabel('$U(r)$ (kcal/mol)',fontsize=size)
ax.legend(loc='upper right',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}_3.png',format='png')
plt.savefig(f'{figname}_3.pdf',format='pdf')
plt.show()