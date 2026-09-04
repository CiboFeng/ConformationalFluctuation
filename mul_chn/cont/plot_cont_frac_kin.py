"""Import Modules"""
from scipy.optimize import curve_fit
from scipy import optimize
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Flx','Mid','Qch']
nmd=len(mdls)
file='cont_frac_kin_annl.pyw'
figname=file[:-4]
nfr=10000

"""Read Data and Calculate"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

t=np.zeros((nmd,nfr))
Qintra=np.zeros((nmd,nfr))
Qinter=np.zeros((nmd,nfr))
for i in range(nmd):
    q=pyw(file,f'({Mdls[i]})')
    t[i]=q[0]
    Qintra[i]=q[1]
    Qinter[i]=q[2]
Qintras=Qintra[:,0:1]
Qintrae=np.mean(Qintra[:,-1000:],axis=-1,keepdims=True)
Qintra=(Qintra-Qintras)/(Qintrae-Qintras)*100
Qinters=Qinter[:,0:1]
Qintere=np.mean(Qinter[:,-1000:],axis=-1,keepdims=True)
Qinter=(Qinter-Qinters)/(Qintere-Qinters)*100

"""Fit Qintra"""
def exp(x,pars):
    y=np.sum(pars[::2])
    n=round(len(pars)/2)
    for i in range(0,2*n,2):
        y+=-pars[i]*np.exp(-pars[i+1]*x)
    return y

nexp_intra=3

init=[1.0,1.0]*nexp_intra
bounds=([-1e1,0.0]*nexp_intra,
        [1e3,1e1]*nexp_intra)
Par_intra=np.ones((nmd,2*nexp_intra))
for i in range(nmd):
    def loss(vars):
        return np.sum((exp(t[i],vars)-Qintra[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par_intra[i]=par
    idx=np.argsort(Par_intra[i,1::2])[::-1]
    Par_intra[i,::2]=Par_intra[i,::2][idx]
    Par_intra[i,1::2]=Par_intra[i,1::2][idx]

Qintra_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    Qintra_fit[i]=exp(t[i],Par_intra[i])
dQintra_fit=Qintra_fit[:,1:]-Qintra_fit[:,:-1]

"""Fit Qinter"""
nexp_inter=2

init=[1.0,1.0]*nexp_inter
bounds=([-1e1,0.0]*nexp_inter,
        [1e3,1e1]*nexp_inter)
Par_inter=np.ones((nmd,2*nexp_inter))
for i in range(nmd):
    def loss(vars):
        return np.sum((exp(t[i],vars)-Qinter[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par_inter[i]=par
    idx=np.argsort(Par_inter[i,1::2])[::-1]
    Par_inter[i,::2]=Par_inter[i,::2][idx]
    Par_inter[i,1::2]=Par_inter[i,1::2][idx]

Qinter_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    Qinter_fit[i]=exp(t[i],Par_inter[i])
dQinter_fit=Qinter_fit[:,1:]-Qinter_fit[:,:-1]

"""Plot"""
fig,axs=plt.subplots(3,2,figsize=(12,12))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.83
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=axs[0,0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i],Qintra[i],color=rgb,linewidth=wth,
            label=(('$-%.0f\\mathrm{e}^{-%.2ft}$'*nexp_intra+'+%.0f')%tuple(list(Par_intra[i])+[np.sum(Par_intra[i,::2])])).replace('--','-'))
    ax.plot(t[i],exp(t[i],Par_intra[i]),color='k',linestyle='--',linewidth=wth)
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
ax.set_xscale('log')
# ax.set_xlim(-4,84)
# ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Contact Formation',fontsize=size)
ax.legend(loc='best',fontsize=0.5*size)

ax=axs[0,1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i],Qinter[i],color=rgb,linewidth=wth,
            label=(('$-%.0f\\mathrm{e}^{-%.2ft}$'*nexp_inter+'+%.0f')%tuple(list(Par_inter[i])+[np.sum(Par_inter[i,::2])])).replace('--','+'))
    ax.plot(t[i],exp(t[i],Par_inter[i]),color='k',linestyle='--',linewidth=wth)
    #         label='%s\n$-%.0f\\mathrm{e}^{-%.2ft}-%.0f\\mathrm{e}^{-%.2ft}+%.0f$'%(mdls[i],Par_inter[i,0],Par_inter[i,3],Par_inter[i,1],Par_inter[i,4],Par_inter[i,2]))
    # ax.plot(t[i],exp_inter(t[i],Par_inter[i,0],Par_inter[i,1],Par_inter[i,2],Par_inter[i,3],Par_inter[i,4]),color='k',linestyle='--',linewidth=wth)
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
ax.set_xscale('log')
# ax.set_yscale('symlog')
# ax.set_ylim(-10,0)
# ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
# ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Contact Formation',fontsize=size)
ax.legend(loc='best',fontsize=0.5*size)

ax=axs[1,0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i,1:1001],dQintra_fit[i,:1000],color=rgb,linewidth=wth)
    ax.plot(t[i,1:],[np.max(dQintra_fit[i,:1000])]*(nfr-1),linewidth=0)
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
ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Rate (ns$^{-1}$)',fontsize=size)

ax=axs[1,1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i,1:1001],dQinter_fit[i,:1000],color=rgb,linewidth=wth)
    ax.plot(t[i,1:],[np.max(dQinter_fit[i,:1000])]*(nfr-1),linewidth=0)
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
ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Rate (ns$^{-1}$)',fontsize=size)

ax=axs[2,0]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i,1:],dQintra_fit[i]/(100-Qintra_fit[i,1:]),color=rgb,linewidth=wth)
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
ax.set_xscale('log')
ax.set_ylim(-0.003,0.06) ### For annl
# ax.set_ylim(-0.005,0.102)
ax.set_xlabel('Intra-chain Time (ns)',fontsize=size)
ax.set_ylabel('Relative Rate (ns$^{-1}$)',fontsize=size)

ax=axs[2,1]
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.plot(t[i,1:],dQinter_fit[i]/(100-Qinter_fit[i,1:]),color=rgb,linewidth=wth)
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
ax.set_xscale('log')
ax.set_ylim(-0.001,0.020) ### For annl
# ax.set_ylim(-0.001,0.018)
ax.set_xlabel('Inter-chain Time (ns)',fontsize=size)
ax.set_ylabel('Relative Rate (ns$^{-1}$)',fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Structure Fluctuation',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

# fig=plt.figure(figsize=(10,5))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# clrs=['r','g','b']

# ax=plt.subplot(1,2,1)
# for i in range(nmd):
#     ax.plot(t[i],Qintra[i],color=clrs[i],linewidth=wth,
#             label=('%s\n'%mdls[i]+('$-%.0f\\mathrm{e}^{-%.2ft}$'*nexp_intra+'+%.0f')%tuple(list(Par_intra[i])+[np.sum(Par_intra[i,::2])])).replace('--','-'))
#     ax.plot(t[i],exp(t[i],Par_intra[i]),color='k',linestyle='--',linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_xlim(-4,84)
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Intra-chain Contact Formation',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

# ax=plt.subplot(1,2,2)
# for i in range(nmd):
#     ax.plot(t[i],Qinter[i],color=clrs[i],linewidth=wth,
#             label=('%s\n'%mdls[i]+('$-%.0f\\mathrm{e}^{-%.2ft}$'*nexp_inter+'+%.0f')%tuple(list(Par_inter[i])+[np.sum(Par_inter[i,::2])])).replace('--','+'))
#     ax.plot(t[i],exp(t[i],Par_inter[i]),color='k',linestyle='--',linewidth=wth)
#     #         label='%s\n$-%.0f\\mathrm{e}^{-%.2ft}-%.0f\\mathrm{e}^{-%.2ft}+%.0f$'%(mdls[i],Par_inter[i,0],Par_inter[i,3],Par_inter[i,1],Par_inter[i,4],Par_inter[i,2]))
#     # ax.plot(t[i],exp_inter(t[i],Par_inter[i,0],Par_inter[i,1],Par_inter[i,2],Par_inter[i,3],Par_inter[i,4]),color='k',linestyle='--',linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_yscale('symlog')
# # ax.set_ylim(-10,0)
# # ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Inter-chain Contact Formation',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

# plt.tight_layout()
# plt.savefig(f'{figname}.png',format='png')
# plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

# """Import Modules"""
# from scipy.optimize import curve_fit
# from scipy import optimize
# import sys
# sys.path.append('D:\\Work\\Code\\Functions')
# from pyw import pyw
# import random
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.pyplot import MultipleLocator
# from matplotlib.collections import LineCollection
# from matplotlib import colors
# from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

# """Set Arguments"""
# mdls=['Flex','Mid','Fix']
# nmd=len(mdls)
# file='cont_frac_kin.pyw'
# figname=file[:-4]
# nfr=10000

# """Read Data and Calculate"""
# t=np.zeros((nmd,nfr))
# Qintra=np.zeros((nmd,nfr))
# Qinter=np.zeros((nmd,nfr))
# for i in range(nmd):
#     q=pyw(file,f'({mdls[i]})')
#     t[i]=q[0]
#     Qintra[i]=q[1]
#     Qinter[i]=q[2]

# """Fit Qintra"""
# nexp_intra=4
# def exp_intra(x,pars):
#     y=pars[-1]
#     n=round((len(pars)-1)/2)
#     for i in range(0,2*n,2):
#         y+=pars[i]*np.exp(-pars[i+1]*x)
#     return y
# init=[1.0,1.0]*nexp_intra+[1.0]
# bounds=([-1e1,0.0]*nexp_intra+[-1e1],
#         [1e3,1e1]*nexp_intra+[1e3])
# Par_intra=np.ones((nmd,2*nexp_intra+1))
# for i in range(nmd):
#     def loss(vars):
#         return np.sum((exp_intra(t[i],vars)-Qintra[i])**2)
#     par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
#     Par_intra[i]=par
#     idx=np.argsort(Par_intra[i,1::2])[::-1]
#     Par_intra[i,::2][:-1]=Par_intra[i,::2][:-1][idx]
#     Par_intra[i,1::2]=Par_intra[i,1::2][idx]

# """Fit Qinter"""
# nexp_inter=4
# def exp_inter(x,pars):
#     y=np.sum(pars[::2])
#     n=round(len(pars)/2)
#     for i in range(0,2*n,2):
#         y+=-pars[i]*np.exp(-pars[i+1]*x)
#     return y
# init=[1.0,1.0]*nexp_inter
# bounds=([-1e1,0.0]*nexp_inter,
#         [1e3,1e1]*nexp_inter)
# Par_inter=np.ones((nmd,2*nexp_inter))
# for i in range(nmd):
#     def loss(vars):
#         return np.sum((exp_inter(t[i],vars)-Qinter[i])**2)
#     par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
#     Par_inter[i]=par
#     idx=np.argsort(Par_inter[i,1::2])[::-1]
#     Par_inter[i,::2]=Par_inter[i,::2][idx]
#     Par_inter[i,1::2]=Par_inter[i,1::2][idx]

# # def exp_inter(x,A,B,C,a,b):
# #     return -A*np.exp(-a*x)-B*np.exp(-b*x)+C
# # init=[1.0,1.0,1.0,1.0,1.0]
# # bounds=([-1e3,-1e3,0.0,0.0,0.0],
# #         [1e3,1e3,50.0,1e1,1e1])
# # # bounds=([-1e3,-1e3,0.0,0.0,0.0],
# # #         [1e3,1e3,50.0,1e1,1e1])
# # Par_inter=np.ones((nmd,5))
# # for i in range(nmd):
# #     par,cov=curve_fit(exp_inter,t[i],Qinter[i],p0=init,bounds=bounds)
# #     # def loss(vars):
# #     #     A,B,C,a,b=vars
# #     #     return np.sum((exp_inter(t[i],A,B,C,a,b)-Qinter[i])**2)
# #     # par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
# #     Par_inter[i]=par
# #     if Par_inter[i,3]<Par_inter[i,4]:
# #         Par_inter[i,:2]=Par_inter[i,:2][::-1]
# #         Par_inter[i,-2:]=Par_inter[i,-2:][::-1]

# """Plot"""
# fig=plt.figure(figsize=(10,5))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# clrs=['r','g','b']

# ax=plt.subplot(1,2,1)
# for i in range(nmd):
#     ax.plot(t[i],Qintra[i],color=clrs[i],linewidth=wth,
#             label=('%s\n'%mdls[i]+('$%.0f\\mathrm{e}^{-%.2ft}$+'*nexp_intra+'%.0f')%tuple(Par_intra[i])).replace('+-','-'))
#     ax.plot(t[i],exp_intra(t[i],Par_intra[i]),color='k',linestyle='--',linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_xlim(-4,84)
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Intra-chain Contact Formation',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

# ax=plt.subplot(1,2,2)
# for i in range(nmd):
#     ax.plot(t[i],Qinter[i],color=clrs[i],linewidth=wth,
#             label=('%s\n'%mdls[i]+('$-%.0f\\mathrm{e}^{-%.2ft}$'*nexp_inter+'+%.0f')%tuple(list(Par_inter[i])+[np.sum(Par_inter[i,::2])])).replace('--','+'))
#     ax.plot(t[i],exp_inter(t[i],Par_inter[i]),color='k',linestyle='--',linewidth=wth)
#     #         label='%s\n$-%.0f\\mathrm{e}^{-%.2ft}-%.0f\\mathrm{e}^{-%.2ft}+%.0f$'%(mdls[i],Par_inter[i,0],Par_inter[i,3],Par_inter[i,1],Par_inter[i,4],Par_inter[i,2]))
#     # ax.plot(t[i],exp_inter(t[i],Par_inter[i,0],Par_inter[i,1],Par_inter[i,2],Par_inter[i,3],Par_inter[i,4]),color='k',linestyle='--',linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.set_xscale('log')
# # ax.set_yscale('symlog')
# # ax.set_ylim(-10,0)
# # ax.set_yticks([-4,-3,-2,-1],['$-10^{4}$','$-10^{3}$','$-10^{2}$','$-10^{1}$'])
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Inter-chain Contact Formation',fontsize=size)
# ax.legend(loc='best',fontsize=0.5*size)

# plt.tight_layout()
# plt.savefig(f'{figname}.png',format='png')
# plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()