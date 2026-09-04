"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
from scipy.optimize import curve_fit
from scipy import optimize
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
file='proc.npz'
figname='chn_clst'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nrep=20
nfr=10000
nch=50

"""Read Data"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

q=np.load(file)
t=q['t']
nclst=q['nclst']
nclst_mean=q['nclst_mean']
nclst_std=q['nclst_std']
treq=q['treq']
treq_mean=q['treq_mean']
treq_std=q['treq_std']
nisol=q['nisol']
nisol_mean=q['nisol_mean']
nisol_std=q['nisol_std']
nsq=q['nsq']
dnsq=q['dnsq']
dnsq_mean=q['dnsq_mean']
dnsq_std=q['dnsq_std']
dnmin=q['dnmin']
dnmin_mean=q['dnmin_mean']
dnmin_std=q['dnmin_std']
ndnmin=np.sum(dnmin!=0,axis=-1)
dtreq=treq[:,:,1:]-treq[:,:,:-1]
dtreq_mean=np.mean(dtreq,axis=1)
# dtreq_mean=treq_mean[:,1:]-treq_mean[:,:-1]

"""Fit nclst"""
def exp1_nclst(x,a):
    return (nch-1)*np.exp(-a*x)+1
Par1_nclst=np.ones(nmd)
for i in range(nmd):
    par,cov=curve_fit(exp1_nclst,t,nclst_mean[i])
    Par1_nclst[i]=par[0]

def exp2_nclst(x,A,a,b):
    return (nch-1)*(A*np.exp(-a*x)+(1-A)*np.exp(-b*x))+1
init=[0.5,1.0,1.0]
bounds=([0.0,0.0,0.0],
        [1.0,np.inf,np.inf])
Par2_nclst=np.ones((nmd,3))
for i in range(nmd):
    par,cov=curve_fit(exp2_nclst,t,nclst_mean[i],p0=init,bounds=bounds)
    Par2_nclst[i]=par
    if Par2_nclst[i,1]<Par2_nclst[i,2]:
        Par2_nclst[i,0]=1-Par2_nclst[i,0]
        Par2_nclst[i,1:]=Par2_nclst[i,1:][::-1]

S0=nch
Sinf=1
P=Par2_nclst[:,0]
k1=Par2_nclst[:,1]
k2=Par2_nclst[:,2]
alf=S0
bet=((P*k1+(1-P)*k2)*Sinf-(1-P)*(k2-k1)*S0)/k1
gma=Sinf
print(bet)

# def exp2_nclst(x,A,B,a,b):
#     return (nch-1)*(A*np.exp(-a*x)+B*np.exp(-b*x))+1
# init=[0.5,0.5,1.0,1.0]
# bounds=([-np.inf,-np.inf,0.0,0.0],
#         [np.inf,np.inf,np.inf,np.inf])
# Par2_nclst=np.ones((nmd,4))
# for i in range(nmd):
#     par,cov=curve_fit(exp2_nclst,t,nclst_mean[i],p0=init,bounds=bounds)
#     Par2_nclst[i]=par

nclst_fit=np.zeros((nmd,nfr))
for i in range(nmd):
    nclst_fit[i]=exp2_nclst(t,Par2_nclst[i,0],Par2_nclst[i,1],Par2_nclst[i,2])
treq_fit=np.zeros((nmd,nch))
for i in range(nmd):
    ks=0
    for j in range(nch):
        for k in range(ks,nfr-1):
            if nclst_fit[i,k+1]<=nch-j<=nclst_fit[i,k]:
                treq_fit[i,nch-j-1]=(nch-j-nclst_fit[i,k])*(t[k+1]-t[k])/(nclst_fit[i,k+1]-nclst_fit[i,k])+t[k]
                ks=k
                break
dtreq_fit=treq_fit[:,1:]-treq_fit[:,:-1]

"""Fit nisol"""
def exp1_nisol(x,a):
    return nch*np.exp(-a*x)
Par1_nisol=np.ones(nmd)
for i in range(nmd):
    par,cov=curve_fit(exp1_nisol,t,nisol_mean[i])
    Par1_nisol[i]=par[0]

def exp2_nisol(x,A,a,b):
    return nch*(A*np.exp(-a*x)+(1-A)*np.exp(-b*x))
init=[0.5,1.0,1.0]
# bounds=([0.0,0.0,0.0],
#         [1.0,np.inf,np.inf])
bounds=([0.0,0.0,0.0],
        [1.0,1e1,1e1])
Par2_nisol=np.ones((nmd,3))
for i in range(nmd):
    # par,cov=curve_fit(exp2_nisol,t,nisol_mean[i],p0=init,bounds=bounds)
    # Par2_nisol[i]=par
    def loss(vars):
        A,a,b=vars
        return np.sum((exp2_nisol(t,A,a,b)-nisol_mean[i])**2)
    par=optimize.differential_evolution(loss,bounds=list(zip(bounds[0],bounds[1]))).x
    Par2_nisol[i]=par
    if Par2_nisol[i,1]<Par2_nisol[i,2]:
        Par2_nisol[i,0]=1-Par2_nisol[i,0]
        Par2_nisol[i,1:]=Par2_nisol[i,1:][::-1]

S0=nch
Sinf=0
P=Par2_nisol[:,0]
k1=Par2_nisol[:,1]
k2=Par2_nisol[:,2]
alf=S0
bet=((P*k1+(1-P)*k2)*Sinf-(1-P)*(k2-k1)*S0)/k1
gma=Sinf
print(bet)

# def exp1(x,a):
#     return np.log(nch*np.exp(-a*np.exp(x)))
# Par=np.ones(nmd)
# for i in range(nmd):
#     par,cov=curve_fit(exp1,np.log(t[1:]),np.log(nclst_mean[i,1:]))
#     Par[i]=par[0]

"""Plot"""
# fig,axs=plt.subplots(2,2,figsize=(12,8))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# cmap_=plt.cm.rainbow
# norm_=colors.Normalize(vmin=0,vmax=nrep)
# color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

# for i in range(nmd):
#     ax=axs[divmod(i,2)]
#     for j in range(nrep):
#         rgb=cmap_(norm_(j))
#         ax.plot(t,nclst[i,j],color=rgb,linewidth=wth,alpha=0.2)
#     ax.plot(t,nclst_mean[i],color='k',linewidth=wth,label='Mean')
#     ax.plot(t,exp1_nclst(t,Par1_nclst[i]),color=(0.7,0.2,0.0),linestyle=':',linewidth=wth,label='$y=49\\mathrm{e}^{-%.2fx}+1$'%Par1_nclst[i])
#     ax.plot(t,exp2_nclst(t,Par2_nclst[i,0],Par2_nclst[i,1],Par2_nclst[i,2]),color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
#             label='$y=49(%.2f\\mathrm{e}^{-%.2fx}+%.2f\\mathrm{e}^{-%.2fx})+1$'%(Par2_nclst[i,0],Par2_nclst[i,1],1-Par2_nclst[i,0],Par2_nclst[i,2]))
#     # ax.plot(t,(nch-1)*(Par2_nclst[i,0]*np.exp(-Par2_nclst[i,2]*t)+Par2_nclst[i,1]*np.exp(-Par2_nclst[i,3]*t))+1,color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
#     #         label='$y=49(%.2f\\mathrm{e}^{-%.2fx}+%.2f\\mathrm{e}^{-%.2fx})+1$'%(Par2_nclst[i,0],Par2_nclst[i,2],Par2_nclst[i,1],Par2_nclst[i,3]))
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xscale('log')
#     # ax.set_yscale('log')
#     ax.set_ylim(-2.5,52.5)
#     if divmod(i,2)[0]==1:
#         ax.set_xlabel('Time (ns)',fontsize=size)
#     ax.set_ylabel('Number of Clusters\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     ax.legend(loc='upper right',fontsize=0.5*size)

# ax=axs[1,1]
# for i in range(nmd):
#     rgb=cmap(i/(nmd-1))
#     ax.plot(t,nclst_mean[i],color=rgb,linewidth=wth)
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
# # ax.set_yscale('log')
# ax.set_ylim(-2.5,52.5)
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Number of Clusters',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Structure Fluctuation',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
# cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_num.png',format='png')
# plt.savefig(f'{figname}_num.pdf',format='pdf')
# plt.show()

# fig,axs=plt.subplots(2,3,figsize=(15,8))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# cmap_=plt.cm.rainbow
# norm_=colors.Normalize(vmin=0,vmax=nrep)
# color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

# for i in range(nmd):
#     ax=axs[0,i]
#     for j in range(nrep):
#         rgb=cmap_(norm_(j))
#         ax.plot(np.arange(1,nch+1),treq[i,j],color=rgb,linewidth=wth,alpha=0.5)
#     ax.plot(np.arange(1,nch+1),treq_mean[i],color='k',linewidth=wth,label='Mean')
#     # ax.plot(t,(nch-1)*np.exp(-Par1[i]*t)+1,color=(0.7,0.2,0.0),linestyle=':',linewidth=wth,label='$y=49\\mathrm{e}^{-%.2fx}+1$'%Par1[i])
#     # ax.plot(t,(nch-1)*(Par2[i,0]*np.exp(-Par2[i,2]*t)+Par2[i,1]*np.exp(-Par2[i,3]*t))/(Par2[i,0]+Par2[i,1])+1,color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
#     #         label='$y=49(%.2f\\mathrm{e}^{-%.2fx}+%.2f\\mathrm{e}^{-%.2fx})+1$'%(Par2[i,0]/(Par2[i,0]+Par2[i,1]),Par2[i,2],Par2[i,1]/(Par2[i,0]+Par2[i,1]),Par2[i,3]))
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     # ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.set_ylim(-2.5,52.5)
#     ax.set_xlabel('Number of Clusters\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     ax.set_ylabel('Required Time (ns)',fontsize=size)
#     if i==0:
#         ax.legend(loc='upper right',fontsize=size)

# ax=axs[1,0]
# for i in range(nmd):
#     rgb=cmap(i/(nmd-1))
#     ax.plot(np.arange(1,nch+1),treq_mean[i],color=rgb,linewidth=wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# # ax.set_xscale('log')
# ax.set_yscale('log')
# # ax.set_ylim(-2.5,52.5)
# ax.set_xlabel('Number of Clusters',fontsize=size)
# ax.set_ylabel('Required Time (ns)',fontsize=size)

# ax=axs[1,1]
# for i in range(nmd):
#     rgb=cmap(i/(nmd-1))
#     ax.errorbar(i,treq_mean[i,24],yerr=treq_std[i,24],
#                 linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
#                 ecolor=rgb,elinewidth=wth,capsize=2*wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# # ax.set_yscale('log')
# ax.set_xticks([0.0,1.0,2.0],np.round(10**rg_flu,2))
# ax.set_xlabel('Structure Fluctuation',fontsize=size)
# ax.set_ylabel('Required Time (ns)\nfor 25 Clusters',fontsize=size)

# ax=axs[1,2]
# for i in range(nmd):
#     rgb=cmap(i/(nmd-1))
#     ax.errorbar(i,treq_mean[i,0],yerr=treq_std[i,0],
#                 linestyle='',marker='.',color=rgb,markersize=5*wth,markeredgewidth=wth,
#                 ecolor=rgb,elinewidth=wth,capsize=2*wth)
# ax.autoscale()
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# # ax.set_yscale('log')
# ax.set_xticks([0.0,1.0,2.0],np.round(10**rg_flu,2))
# ax.set_xlabel('Structure Fluctuation',fontsize=size)
# ax.set_ylabel('Required Time (ns)\nfor 1 Clusters',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Structure Fluctuation',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
# cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_time.png',format='png')
# plt.savefig(f'{figname}_time.pdf',format='pdf')
# plt.show()

fig,axs=plt.subplots(1,1,figsize=(8,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.8
cbarwth=0.03
cmap_=plt.cm.rainbow
norm_=colors.Normalize(vmin=0,vmax=nrep)
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=axs
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    # ax.plot(treq_mean[i,1:],-1/dtreq_mean[i],color=rgb,linewidth=wth)
    ax.plot(treq_fit[i,1:],-1/dtreq_fit[i],color=rgb,linewidth=wth)
    # ax.plot(treq_fit[i],color=rgb,linewidth=wth)
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
# ax.set_ylim(-2.5,52.5)
ax.set_xlabel('Time (ns)',fontsize=size)
ax.set_ylabel('Rate (ns$^{-1}$)',fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs.get_position().y1
bot=axs.get_position().y0
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

plt.savefig(f'{figname}_rate.png',format='png')
plt.savefig(f'{figname}_rate.pdf',format='pdf')
plt.show()

# fig,axs=plt.subplots(2,2,figsize=(12,8))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# cmap_=plt.cm.rainbow
# norm_=colors.Normalize(vmin=0,vmax=nrep)
# color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
# nodes=[0/3,1/3,2/3,3/3]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

# for i in range(nmd):
#     ax=axs[divmod(i,2)]
#     for j in range(nrep):
#         rgb=cmap_(norm_(j))
#         ax.plot(t,nisol[i,j],color=rgb,linewidth=wth,alpha=0.2)
#     ax.plot(t,nisol_mean[i],color='k',linewidth=wth,label='Mean')
#     ax.plot(t,exp1_nisol(t,Par1_nisol[i]),color=(0.7,0.2,0.0),linestyle=':',linewidth=wth,label='$y=50\\mathrm{e}^{-%.2fx}$'%Par1_nisol[i])
#     ax.plot(t,exp2_nisol(t,Par2_nisol[i,0],Par2_nisol[i,1],Par2_nisol[i,2]),color=(0.7,0.2,0.0),linestyle='--',linewidth=wth,
#             label='$y=50(%.2f\\mathrm{e}^{-%.2fx}+%.2f\\mathrm{e}^{-%.2fx})$'%(Par2_nisol[i,0],Par2_nisol[i,1],1-Par2_nisol[i,0],Par2_nisol[i,2]))
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xscale('log')
#     # ax.set_yscale('log')
#     ax.set_ylim(-2.5,52.5)
#     if divmod(i,2)[0]==1:
#         ax.set_xlabel('Time (ns)',fontsize=size)
#     ax.set_ylabel('Number of Isolated\nChains of Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     ax.legend(loc='upper right',fontsize=0.5*size)

# ax=axs[1,1]
# for i in range(nmd):
#     rgb=cmap(i/(nmd-1))
#     ax.plot(t,nisol_mean[i],color=rgb,linewidth=wth)
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
# # ax.set_yscale('log')
# ax.set_ylim(-2.5,52.5)
# ax.set_xlabel('Time (ns)',fontsize=size)
# ax.set_ylabel('Number of Isolated Chains',fontsize=size)

# fig.subplots_adjust(right=figwth)
# plt.tight_layout(rect=[0,0,figwth,1])
# top=axs[0,-1].get_position().y1
# bot=axs[-1,-1].get_position().y0
# cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
# cmap_disc=ListedColormap(cmap(np.linspace(0,1,nmd)))
# sm=plt.cm.ScalarMappable(cmap=cmap_disc)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.set_label('Structure Fluctuation',fontsize=size)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
# cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
# cbar.set_ticks(np.linspace(1/(2*nmd),1-1/(2*nmd),nmd))
# cbar.set_ticklabels([f'{10**rg_flu[i]:.2f}' for i in range(nmd)])
# cbar.outline.set_linewidth(wth)

# plt.savefig(f'{figname}_isol.png',format='png')
# plt.savefig(f'{figname}_isol.pdf',format='pdf')
# plt.show()

# fig=plt.figure(figsize=(15,12))
# wth=2
# size=20
# lenmaj=15
# lenmin=8
# xtick=0.1
# ytick=20
# figwth=0.87
# cbarwth=0.03
# # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# # nodes=[0/4,1/4,2/4,3/4,4/4]
# # cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# cmap=plt.cm.rainbow
# norm=colors.Normalize(vmin=0,vmax=nrep)
# lbls=['Flex','Mid','Fix']

# for i in range(nmd):
#     ax=plt.subplot(3,3,i+1)
#     for j in range(nrep):
#         rgb=cmap(norm(j))
#         ax.plot(t,np.sqrt(nsq[i,j]),color=rgb,linewidth=wth,alpha=0.5)
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xscale('log')
#     # ax.set_yscale('log')
#     ax.set_ylim(4.8,52.2)
#     # ax.set_xlabel('Time (ns) of Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     if i==0:
#         ax.set_ylabel('Clusters Amount',fontsize=size)

# for i in range(nmd):
#     ax=plt.subplot(3,3,nmd+i+1)
#     for j in range(nrep):
#         rgb=cmap(norm(j))
#         ax.scatter(t[:-1],dnsq[i,j],color=rgb,s=10*wth,alpha=0.5)
#     # ax.fill_between(t[:-1],dnsq_mean[i]-dnsq_std[i],dnsq_mean[i]+dnsq_std[i],color='k',alpha=0.5)
#     # ax.plot(t[1:],dnsq_mean[i],color='k',linewidth=wth)
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xscale('log')
#     # ax.set_yscale('log')
#     ax.set_ylim(-1.3,28.3)
#     ax.set_xlabel('Time (ns) of Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     if i==0:
#         ax.set_ylabel('Merging Amount',fontsize=size)

# for i in range(nmd):
#     ax=plt.subplot(3,3,2*nmd+i+1)
#     for j in range(nrep):
#         rgb=cmap(norm(j))
#         for k in range(np.max(ndnmin[i,j,:])):
#             ax.scatter(t[:-1],dnmin[i,j,:,k],color=rgb,s=10*wth,alpha=0.5)
#         ax.plot(t[:-1],ndnmin[i,j,:],color=rgb,linewidth=wth,alpha=0.5)
#     # ax.fill_between(t[:-1],dnsq_mean[i]-dnsq_std[i],dnsq_mean[i]+dnsq_std[i],color='k',alpha=0.5)
#     # ax.plot(t[1:],dnsq_mean[i],color='k',linewidth=wth)
#     ax.autoscale()
#     ax.minorticks_on()
#     ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.spines['bottom'].set_linewidth(wth)
#     ax.spines['top'].set_linewidth(wth)
#     ax.spines['left'].set_linewidth(wth)
#     ax.spines['right'].set_linewidth(wth)
#     ax.set_xscale('log')
#     # ax.set_yscale('log')
#     ax.set_ylim(-1.3,26.3)
#     ax.set_xlabel('Time (ns) of Flu.=%.2f'%10**rg_flu[i],fontsize=size)
#     if i==0:
#         ax.set_ylabel('Merging Amount Minimum',fontsize=size)

# plt.tight_layout()
# plt.savefig(f'{figname}_amt.png',format='png')
# plt.savefig(f'{figname}_amt.pdf',format='pdf')
# plt.show()