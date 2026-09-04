"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
dir0='../../max_entr/rg/rg_heat'
dir='dst'
figname1=dir
figname2='fs_dgrm'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=1000
unit=0.0738e-3 ### 1g/(mol*AA)=0.0738mg/ml
nintp=10000
bet=0.325

"""Read Data and Calculate"""
rg_mean=np.zeros(nmd)
rg_std=np.zeros(nmd)
for i in range(nmd):
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Mean')
    rg_mean[i]=q[0,0]
    q=pyw(f'{dir0}/{mdls[i]}_300.0.pyw','Deviation')
    rg_std[i]=q[0,0]
rg_flu=np.round(np.log10(rg_std/rg_mean),1)

z=np.zeros((nmd,nT,nhi))
dst=np.zeros((nmd,nT,nhi))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir}/{mdls[i]}_{Ts[j]}.pyw','Density:')
        z[i,j]=q[0]
        dst[i,j]=q[1]

def rho(z,z1,z2,d,rhol,rhoh):
    return rhol+(rhoh-rhol)/2*(np.tanh((z-z1)/d)-np.tanh((z-z2)/d))
init1=[np.min(z)+0.3*(np.max(z)-np.min(z)),np.min(z)+0.7*(np.max(z)-np.min(z)),1.0,np.min(dst),np.max(dst)]
bounds1=([np.min(z),np.min(z),0.0,0.0,0.0],
        [np.max(z),np.max(z),np.inf,np.inf,np.inf])

def rho_dlt(x,A,Tc):
    return A*(Tc-x)**bet
init2=[1.0,np.max(Ts)*2]
bounds2=([0.0,np.max(Ts)],
        [np.inf,np.inf])

def rho_sum(x,B,rhoc):
    return B*x+2*rhoc

def rholh(x,A,B,Tc,rhoc):
    n=len(x)//2
    x1=x[:n]
    x2=x[n:]
    y1=rhoc-A/2*(Tc-x1)**bet+B/2*(Tc-x1)
    y2=rhoc+A/2*(Tc-x2)**bet+B/2*(Tc-x2)
    return np.concatenate([y1,y2])
init4=[1.0,1.0,np.max(Ts)*2,1.0]
bounds4=([0.0,0.0,np.max(Ts),0.0],
        [np.inf,np.inf,np.inf,np.inf])

rhol=np.zeros((nmd,nT))
rhoh=np.zeros((nmd,nT))
A=np.ones(nmd)
Tc=np.ones(nmd)
B=np.ones(nmd)
rhoc=np.ones(nmd)
Tintp=np.zeros((nmd,nintp))
rholintp=np.zeros((nmd,nintp))
rhohintp=np.zeros((nmd,nintp))
for i in range(nmd):
    for j in range(nT):
        par,cov=curve_fit(rho,z[i,j],dst[i,j],p0=init1,bounds=bounds1) ### Not only one set optimal parameters
        rhol[i,j],rhoh[i,j]=par[-2:]
    # par,cov=curve_fit(rho_dlt,Ts,rhoh[i]-rhol[i],p0=init2,bounds=bounds2)
    # A[i],Tc[i]=par
    # par,cov=curve_fit(rho_sum,par[1]-Ts,rhoh[i]+rhol[i])
    # B[i],rhoc[i]=par
    par,cov=curve_fit(rholh,np.concatenate([Ts,Ts]),np.concatenate([rhol[i],rhoh[i]]),p0=init4,bounds=bounds4)
    A[i],B[i],Tc[i],rhoc[i]=par
    Tintp[i]=np.logspace(np.log10(Tc[i]),np.log10(np.min(Ts)),nintp)[::-1]
    rholintp[i]=rhoc[i]-A[i]/2*(Tc[i]-Tintp[i])**bet+B[i]/2*(Tc[i]-Tintp[i])
    rhohintp[i]=rhoc[i]+A[i]/2*(Tc[i]-Tintp[i])**bet+B[i]/2*(Tc[i]-Tintp[i])
    for j in range(nintp):
        if Tintp[i,j]<Ts[-1]:
            Tleft=max(x for x in Ts if x<=Tintp[i,j])
            Tright=min(x for x in Ts if x>=Tintp[i,j])
            rholeft=rhol[i,Ts.index(Tleft)]
            rhoright=rhol[i,Ts.index(Tright)]
            rholintp[i,j]=(Tintp[i,j]-Tleft)*(rhoright-rholeft)/(Tright-Tleft)+rholeft
            rholeft=rhoh[i,Ts.index(Tleft)]
            rhoright=rhoh[i,Ts.index(Tright)]
            rhohintp[i,j]=(Tintp[i,j]-Tleft)*(rhoright-rholeft)/(Tright-Tleft)+rholeft

"""Plot"""
fig,axs=plt.subplots(nT,2,figsize=(12,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.86
cbarwth=0.03
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

for i in range(nT):
    ax=axs[i,0]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(z[j,i],dst[j,i]*unit,color=rgb,linewidth=wth)
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
    ax.set_ylim(-0.2,2.2)
    if i==nT-1:
        ax.set_xlabel('Z ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Density (g/ml)\nat %i K'%Ts[i],fontsize=size)
    
    ax=axs[i,1]
    for j in range(nmd):
        rgb=cmap(j/(nmd-1))
        ax.plot(z[j,i],dst[j,i]*unit,color=rgb,linewidth=wth)
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
    ax.set_yscale('log')
    ax.set_ylim(5e-7,5)
    if i==nT-1:
        ax.set_xlabel('Z ($\\mathrm{\\AA}$)',fontsize=size)

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

plt.savefig(f'{figname1}_1.png',format='png')
plt.savefig(f'{figname1}_1.pdf',format='pdf')
plt.show()

fig,axs=plt.subplots(nmd,2,figsize=(12,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
figwth=0.86
cbarwth=0.03
color=['#087F5B','#0CA678','#20C997','#63E6BE']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

for i in range(nmd):
    ax=axs[i,0]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(z[i,j],dst[i,j]*unit,color=rgb,linewidth=wth)
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
    ax.set_ylim(-0.2,2.2)
    if i==nmd-1:
        ax.set_xlabel('Z ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Density (g/ml)\nof Flu.=%.2f'%10**rg_flu[i],fontsize=size)
    
    ax=axs[i,1]
    for j in range(nT):
        rgb=cmap(j/(nT-1))
        ax.plot(z[i,j],dst[i,j]*unit,color=rgb,linewidth=wth)
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
    ax.set_yscale('log')
    ax.set_ylim(5e-7,5)
    if i==nmd-1:
        ax.set_xlabel('Z ($\\mathrm{\\AA}$)',fontsize=size)

fig.subplots_adjust(right=figwth)
plt.tight_layout(rect=[0,0,figwth,1])
top=axs[0,-1].get_position().y1
bot=axs[-1,-1].get_position().y0
cbar_ax=fig.add_axes([figwth+cbarwth,bot,cbarwth,top-bot])
cmap_disc=ListedColormap(cmap(np.linspace(0,1,nT)))
sm=plt.cm.ScalarMappable(cmap=cmap_disc)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Temperature (K)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',labelsize=size,direction='in',width=wth,length=lenmin)
cbar.ax.tick_params(which='minor',labelsize=size,direction='in',width=wth,length=0)
cbar.set_ticks(np.linspace(1/(2*nT),1-1/(2*nT),nT))
cbar.set_ticklabels([int(Ts[i]) for i in range(nT)])
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname1}_2.png',format='png')
plt.savefig(f'{figname1}_2.pdf',format='pdf')
plt.show()

fig=plt.figure(figsize=(6,6))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
color=['#105186','#3D86B6','#67B3DA','#9DDAF2']
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))

ax=plt.subplot(111)
for i in range(nmd):
    rgb=cmap(i/(nmd-1))
    ax.scatter(rhol[i]*unit,Ts,color=rgb,linewidth=wth,label='Flu.=%.2f'%10**rg_flu[i],zorder=i*4)
    ax.scatter(rhoh[i]*unit,Ts,color=rgb,linewidth=wth,zorder=i*4+1)
    ax.plot(rholintp[i]*unit,Tintp[i],color=rgb,linewidth=wth,zorder=i*4+2)
    ax.plot(rhohintp[i]*unit,Tintp[i],color=rgb,linewidth=wth,zorder=i*4+3)
ax.scatter(rhoc*unit,Tc,c='k',marker='*',linewidth=wth,zorder=12)
# for i in range(nmd):
#     ax.scatter(Tc[i]-Ts,rhoh[i]-rhol[i],color=rgb,linewidth=wth,label='Flu.=%.2f'%10**rg_flu[i])
#     ax.scatter(Tc[i]-Ts,rhoh[i]+rhol[i],color=rgb,linewidth=wth)
#     ax.plot(Tc[i]-Tintp[i],rhohintp[i]-rholintp[i],color=rgb,linewidth=wth)
#     ax.plot(Tc[i]-Tintp[i],rhohintp[i]+rholintp[i],color=rgb,linewidth=wth)
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
ax.set_xlabel('Density (g/ml)',fontsize=size)
ax.set_ylabel('Temperature (K)',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.legend(loc='best',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname2}.png',format='png')
plt.savefig(f'{figname2}.pdf',format='pdf')
plt.show()