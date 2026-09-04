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
dir1='../../dst/dst'
file2='../../strs/surf_tens.npz'
dir2='../../rg/rg_slab'
dir3='../../../max_entr/rg/rg_heat'
figname='thermo'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=1000
unit1=0.0738e-3 ### 1g/(mol*AA)=0.0738mg/ml
unit2=0.101325 ### mN/m
nintp=10000
bet=0.325
nhiz=250
nhirg=50
nhirg0=200
nsms=10

"""Density"""
z=np.zeros((nmd,nT,nhi))
dst=np.zeros((nmd,nT,nhi))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Density:')
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

"""Surface Tension"""
q=np.load(file2)
F=q['F']
p_F=q['p_F']
F_mean=q['F_mean']
F_std=q['F_std']

"""Radius of Gyration"""
z_rg=np.zeros((nmd,nT,nhiz))
rg=np.zeros((nmd,nT,nhirg))
p=np.zeros((nmd,nT,nhiz,nhirg))
rg_mean=np.zeros((nmd,nT,nhiz))
rg_std=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir2}/{mdls[i]}_{Ts[j]}.npz')
        z_rg[i,j]=q['z']
        rg[i,j]=q['rg']
        p[i,j]=q['p']
        rg_mean[i,j]=q['rg_mean']
        rg_std[i,j]=q['rg_std']

pz=np.sum(p,axis=-1)
pz/=np.sum(pz*(np.max(z_rg)-np.min(z_rg))/nhiz,axis=-1,keepdims=True)
z1=np.zeros((nmd,nT))
z2=np.zeros((nmd,nT))
d=np.zeros((nmd,nT))
pedge=np.zeros((nmd,nT,3,nhirg))
rgedge_mean=np.zeros((nmd,nT,3))
rgedge_std=np.zeros((nmd,nT,3))
for i in range(nmd):
    for j in range(nT):
        z1[i,j],z2[i,j],d[i,j]=pyw(f'{dir1}/{mdls[i]}_{Ts[j]}.pyw','Position')[:3,0]
        idxz=[[],[],[]]
        idxz[0]+=[round(nhiz/2)-1,round(nhiz/2)]
        for k in range(nhiz):
            if z1[i,j]+d[i,j]<z_rg[i,j,k]<z2[i,j]-d[i,j]:
                idxz[1]+=[k]
            if z1[i,j]-d[i,j]<z_rg[i,j,k]<z1[i,j]+d[i,j] or z2[i,j]-d[i,j]<z[i,j,k]<z2[i,j]+d[i,j]:
                idxz[2]+=[k]
        for k in range(len(idxz)):
            idx=(i,j,idxz[k])
            pedge[i,j,k]=np.sum(p[idx],axis=0)
            rgedge_mean[i,j,k]=np.sum(pz[idx]*rg_mean[idx])/np.sum(pz[idx])
            rgedge_std[i,j,k]=np.sqrt(np.sum(pz[idx]*rg_std[idx]**2)/np.sum(pz[idx])+ \
                    np.sum(pz[idx]*(rg_mean[idx]-rgedge_mean[i,j,k])**2)/np.sum(pz[idx]))
pedge/=np.sum(pedge*(np.max(rg)-np.min(rg))/nhirg,axis=-1,keepdims=True)

rg_dev=np.zeros((nmd,nT,3))
rg_dev[:,:,0]=rgedge_mean[:,:,1]/np.mean(rgedge_mean[:,:,1])-1
rg_dev[:,:,1]=rgedge_std[:,:,1]/np.mean(rgedge_std[:,:,1])-1
rg_dev[:,:,2]=rgedge_std[:,:,1]/rgedge_mean[:,:,1]/np.mean(rgedge_std[:,:,1]/rgedge_mean[:,:,1])-1

rg_mean[rg_mean==0.0]=np.nan
rg_std[rg_std==0.0]=np.nan

rg_mean_sms=np.zeros((nmd,nT,nhiz))
rg_std_sms=np.zeros((nmd,nT,nhiz))
for i in range(nmd):
    for j in range(nT):
        for k in range(nhiz):
            idx=(i,j,list(np.arange(max(0,k-round(nsms/2)),min(nhiz,k+round(nsms/2)))))
            rg_mean_sms[i,j,k]=np.nansum(pz[idx]*rg_mean[idx])/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0))
            rg_std_sms[i,j,k]=np.sqrt(np.nansum(pz[idx]*rg_std[idx]**2)/np.nansum(pz[idx]*np.heaviside(rg_std[idx],0.0))+ \
                    np.nansum(pz[idx]*(rg_mean[idx]-rg_mean_sms[i,j,k])**2)/np.nansum(pz[idx]*np.heaviside(rg_mean[idx],0.0)))

rg0=np.zeros((nmd,nT,nhirg0))
p_rg0=np.zeros((nmd,nT,nhirg0))
rg0_mean=np.zeros((nmd,nT))
rg0_std=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        q=pyw(f'{dir3}/{mdls[i]}_{Ts[j]}.pyw','Probability')
        rg0[i,j]=q[0]
        p_rg0[i,j]=q[1]
        q=pyw(f'{dir3}/{mdls[i]}_{Ts[j]}.pyw','Mean')
        rg0_mean[i,j]=q[0,0]
        q=pyw(f'{dir3}/{mdls[i]}_{Ts[j]}.pyw','Deviation')
        rg0_std[i,j]=q[0,0]

"""Plot"""
fig=plt.figure(figsize=(25,15))
gs=fig.add_gridspec(nrows=1,ncols=3,width_ratios=[2,1,1],wspace=0.3)
Gs=[gs[0].subgridspec(nrows=3,ncols=1,height_ratios=[1,1,1],hspace=0.5),
    gs[1].subgridspec(nrows=2,ncols=1,height_ratios=[1,1],hspace=0.3),
    gs[2].subgridspec(nrows=4,ncols=1,height_ratios=[1,1,1,1],hspace=0.5)]
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
ann=[['(A)','(B)','(C)'],['(D)','(E)'],['(F)','(G)','(H)','(I)']]
locs=['Center','Inside','Surface','Isolated']
qs=['$\\langle{}R_\\mathrm{g}\\rangle$','std.$(R_\\mathrm{g})$','cv.$(R_\\mathrm{g})$']

for i in range(nmd):
    ax=axs[0][i]
    ax.annotate(ann[0][i],xy=(0.0,1.075),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel(Mdls[i],fontdict=font)

ax=axs[1][0]
ax.annotate(ann[1][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.scatter(rhol[i]*unit1,Ts,color=clrs[i],linewidth=wth,zorder=i*4)
    ax.scatter(rhoh[i]*unit1,Ts,color=clrs[i],linewidth=wth,zorder=i*4+1)
    ax.plot(rholintp[i]*unit1,Tintp[i],color=clrs[i],linewidth=wth,label=Mdls[i],zorder=i*4+2)
    ax.plot(rhohintp[i]*unit1,Tintp[i],color=clrs[i],linewidth=wth,zorder=i*4+3)
ax.scatter(rhoc*unit1,Tc,c='k',marker='*',s=100*wth,label='$T_\\mathrm{c}$',zorder=12)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xlabel('$\\rho$ (g/ml)',fontdict=font)
ax.set_ylabel('$T$ (K)',fontdict=font)
leg=ax.legend(loc='best',handlelength=0.0,prop={'family':font['family'],'size':font['size']})
for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
    if hasattr(hnd,'get_color'):
        legc=hnd.get_color()
    elif hasattr(hnd,'get_facecolor'):
        legc=hnd.get_facecolor()
    else:
        legc='k' 
    txt.set_color(legc)

ax=axs[1][1]
ax.annotate(ann[1][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
spc=3
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,F_mean[i,:]*unit2,yerr=F_std[i,:]*unit2/10,
                     linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                     ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nT)*(nmd+spc)+(nmd-1)/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$\\Gamma$ (mN/m)',fontdict=font)

spc=3
for i in range(1,3):
    ax=axs[2][i-1]
    ax.annotate(ann[2][i-1],xy=(0.0,1.1),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    for j in range(nmd):
        ax.errorbar(np.arange(nT)*(nmd+spc)+j,rgedge_mean[j,:,i],yerr=rgedge_std[j,:,i],
                        linestyle='',marker='.',color=clrs[j],markersize=5*wth,markeredgewidth=wth,
                        ecolor=clrs[j],elinewidth=wth,capsize=2*wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(np.arange(nT)*(nmd+spc)+nmd/2,np.array(Ts).astype(int),fontdict=font)
    ax.set_xlabel('$T$ (K)',fontdict=font)
    ax.set_ylabel('$R\\mathrm{_g~(\\AA)}$\nat %s'%locs[i],fontdict=font)

ax=axs[2][2]
ax.annotate(ann[2][2],xy=(0.0,1.1),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
for i in range(nmd):
    ax.errorbar(np.arange(nT)*(nmd+spc)+i,rg0_mean[i,:],yerr=rg0_std[i,:],
                    linestyle='',marker='.',color=clrs[i],markersize=5*wth,markeredgewidth=wth,
                    ecolor=clrs[i],elinewidth=wth,capsize=2*wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ax.set_xticks(np.arange(nT)*(nmd+spc)+nmd/2,np.array(Ts).astype(int),fontdict=font)
ax.set_xlabel('$T$ (K)',fontdict=font)
ax.set_ylabel('$R\\mathrm{_g~(\\AA)}$',fontdict=font)

ylims=[axs[2][i].get_ylim() for i in range(3)]
for i in range(3):
    axs[2][i].set_ylim(np.min(ylims),np.max(ylims))

ax=axs[2][3]
ax.annotate(ann[2][3],xy=(0.0,1.1),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
ymin=1.1*np.min(rg_dev)-0.1*np.max(rg_dev)
for i in range(nmd):
    ax.bar(np.arange(3)*(nmd+spc)+i+1,rg_dev[i,0]-ymin,width=1.0,color=clrs[i],edgecolor='k',linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
ax.xaxis.set_minor_locator(AutoMinorLocator(0))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_yscale('log')
ax.set_xticks(np.arange(3)*(nmd+spc)+nmd/2,qs,rotation=30,fontdict=font)
yticks=np.array([-0.1,0.0,0.1])
ax.set_yticks(yticks-ymin,yticks,fontdict=font)
ax.set_ylabel('Relative Deviation',fontdict=font)

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