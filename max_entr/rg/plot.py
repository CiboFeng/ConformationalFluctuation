"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
file0='../../all_atom/rg/rg.pyw'
file1='../../hps/rg/rg.pyw'

csgm=[1.5,1.55,1.6]
esgm=1.0
itr=[20,20,20,20]
nc=len(csgm)

# csgm=[1.37,1.4,1.5]
# esgm=0.6
# itr=[20,20,20]
# nc=len(csgm)

# csgm=1.5
# esgm=[0.6,0.8,0.82,1.0]
# itr=[20,20,20,20]
# ne=len(esgm)

file='fus_sgm/rg.pyw'
figname=file.split('/')[-1][:-4]+'_sgm'
nhi=200

"""Read Data and Calculate"""
q=pyw(file0,'Probability')
rg0=q[0]
p_rg0=q[1]
q=pyw(file0,'Mean')
rg0_mean=q[0][0]
q=pyw(file0,'Deviation')
rg0_std=q[0][0]

q=pyw(file1,'Probability')
rg1=q[0]
p_rg1=q[1]
q=pyw(file1,'Mean')
rg1_mean=q[0][0]
q=pyw(file1,'Deviation')
rg1_std=q[0][0]

if type(csgm)==list:
    rg=np.zeros((nc,nhi))
    p_rg=np.zeros((nc,nhi))
    rg_mean=np.zeros(nc)
    rg_std=np.zeros(nc)
    for i in range(nc):
        q=pyw(file.replace('sgm',f'{csgm[i]}sgm_{esgm}sgm').replace('rg',f'rg_itr{itr[i]}'),'Probability')
        rg[i]=q[0]
        p_rg[i]=q[1]
        q=pyw(file.replace('sgm',f'{csgm[i]}sgm_{esgm}sgm').replace('rg',f'rg_itr{itr[i]}'),'Mean')
        rg_mean[i]=q[0][0]
        q=pyw(file.replace('sgm',f'{csgm[i]}sgm_{esgm}sgm').replace('rg',f'rg_itr{itr[i]}'),'Deviation')
        rg_std[i]=q[0][0]
else:
    rg=np.zeros((ne,nhi))
    p_rg=np.zeros((ne,nhi))
    rg_mean=np.zeros(ne)
    rg_std=np.zeros(ne)
    for i in range(ne):
        q=pyw(file.replace('sgm',f'{csgm}sgm_{esgm[i]}sgm').replace('rg',f'rg_itr{itr[i]}'),'Probability')
        rg[i]=q[0]
        p_rg[i]=q[1]
        q=pyw(file.replace('sgm',f'{csgm}sgm_{esgm[i]}sgm').replace('rg',f'rg_itr{itr[i]}'),'Mean')
        rg_mean[i]=q[0][0]
        q=pyw(file.replace('sgm',f'{csgm}sgm_{esgm[i]}sgm').replace('rg',f'rg_itr{itr[i]}'),'Deviation')
        rg_std[i]=q[0][0]

"""Plot"""
if type(csgm)==list:
    wth=2
    size=20
    lenmaj=15
    lenmin=8
    figrig=0.82
    cbarwth=0.03
    xtick=0.1
    ytick=20
    # ytick=0.002
    # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
    # nodes=[0/4,1/4,2/4,3/4,1/4]
    # cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    cmap=plt.cm.rainbow
    norm=colors.Normalize(vmin=min(csgm),vmax=max(csgm))

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    ax=axs[0]
    for i in range(nc):
        normc=norm(csgm[i])
        rgb=cmap(normc)
        ax.plot(rg[i],p_rg[i],color=rgb,linewidth=wth)
    ax.plot(rg0,p_rg0,'k',linewidth=wth,label='AA')
    ax.plot(rg1,p_rg1,'k:',linewidth=wth,label='HPS')
    # ax.plot([rg0_mean,rg0_mean],[np.min(p_rg0),np.max(p_rg0)],'--',linewidth=wth,label='Mean')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(10,28)
    ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Probability Density ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
    ax.legend(loc='best',fontsize=size)

    ax=axs[1]
    ax.plot(csgm,rg_mean,'--',color=(0.5,0.5,0.5),linewidth=wth)
    ax.fill_between([np.min(csgm),np.max(csgm)],[rg0_mean-rg0_std,rg0_mean-rg0_std],[rg0_mean+rg0_std,rg0_mean+rg0_std],color='k',linewidth=0,alpha=0.2)
    ax.plot([np.min(csgm),np.max(csgm)],[rg0_mean,rg0_mean],'k',linewidth=wth)
    ax.fill_between([np.min(csgm),np.max(csgm)],[rg1_mean-rg1_std,rg1_mean-rg1_std],[rg1_mean+rg1_std,rg1_mean+rg1_std],color='k',linewidth=0,alpha=0.2)
    ax.plot([np.min(csgm),np.max(csgm)],[rg1_mean,rg1_mean],'k:',linewidth=wth)
    for i in range(nc):
        normc=norm(csgm[i])
        rgb=cmap(normc)
        ax.errorbar(csgm[i],rg_mean[i],yerr=rg_std[i],
                     marker='.',markersize=0,markeredgewidth=2*wth,
                     ecolor=(0.5,0.5,0.5),elinewidth=2*wth,capsize=2.5*wth)
        ax.errorbar(csgm[i],rg_mean[i],yerr=rg_std[i],
                     marker='.',markersize=0,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
        ax.plot(csgm[i],rg_mean[i],'.',color=(0.5,0.5,0.5),markersize=8*wth)
        ax.plot(csgm[i],rg_mean[i],'.',color=rgb,markersize=5*wth)
    ax.plot([],[],label='Exclusion Radius\n= %s $\\sigma_{ij}$'%esgm)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('Capture Radius ($\\sigma_{ij}$)',fontsize=size)
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
    ax.legend(loc='center',fontsize=size,handlelength=0.0)

    fig.subplots_adjust(right=figrig)
    plt.tight_layout(rect=[0,0,figrig,1])
    top=ax[-1].get_position().y1
    bot=ax[-1].get_position().y0
    cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
    sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm,cax=cbar_ax)
    cbar.set_label('Capture Radius ($\\sigma_{ij}$)',fontsize=size)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(0.05))
    # cbar.set_ticks(csgm)
    cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
    cbar.outline.set_linewidth(wth)

    plt.savefig(f'{figname}_{esgm}sgm.png',format='png')
    plt.savefig(f'{figname}_{esgm}sgm.pdf',format='pdf')
    plt.show()

else:
    wth=2
    size=20
    lenmaj=15
    lenmin=8
    xtick=0.1
    ytick=20
    figrig=0.82
    cbarwth=0.03
    # ytick=0.002
    # color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
    # nodes=[0/4,1/4,2/4,3/4,4/4]
    # cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    cmap=plt.cm.rainbow
    norm=colors.Normalize(vmin=min(esgm),vmax=max(esgm))

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    ax=axs[0]
    for i in range(ne):
        normc=norm(esgm[i])
        rgb=cmap(normc)
        ax.plot(rg[i],p_rg[i],color=rgb,linewidth=wth)
    ax.plot(rg0,p_rg0,'k',linewidth=wth,label='AA')
    ax.plot(rg1,p_rg1,'k:',linewidth=wth,label='HPS')
    # ax.plot([rg0_mean,rg0_mean],[np.min(p_rg0),np.max(p_rg0)],'--',linewidth=wth,label='Mean')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlim(10,28)
    ax.set_xlabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
    ax.set_ylabel('Probability Density ($\\mathrm{\\AA}^{-1}$)',fontsize=size)
    ax.legend(loc='best',fontsize=size)

    ax=axs[1]
    ax.plot(esgm,rg_mean,'--',color=(0.5,0.5,0.5),linewidth=wth)
    ax.fill_between([np.min(esgm),np.max(esgm)],[rg0_mean-rg0_std,rg0_mean-rg0_std],[rg0_mean+rg0_std,rg0_mean+rg0_std],color='k',linewidth=0,alpha=0.2)
    ax.plot([np.min(esgm),np.max(esgm)],[rg0_mean,rg0_mean],'k',linewidth=wth)
    ax.fill_between([np.min(esgm),np.max(esgm)],[rg1_mean-rg1_std,rg1_mean-rg1_std],[rg1_mean+rg1_std,rg1_mean+rg1_std],color='k',linewidth=0,alpha=0.2)
    ax.plot([np.min(esgm),np.max(esgm)],[rg1_mean,rg1_mean],'k:',linewidth=wth)
    for i in range(ne):
        normc=norm(esgm[i])
        rgb=cmap(normc)
        ax.errorbar(esgm[i],rg_mean[i],yerr=rg_std[i],
                     marker='.',markersize=0,markeredgewidth=2*wth,
                     ecolor=(0.5,0.5,0.5),elinewidth=2*wth,capsize=2.5*wth)
        ax.errorbar(esgm[i],rg_mean[i],yerr=rg_std[i],
                     marker='.',markersize=0,markeredgewidth=wth,
                     ecolor=rgb,elinewidth=wth,capsize=2*wth)
        ax.plot(esgm[i],rg_mean[i],'.',color=(0.5,0.5,0.5),markersize=8*wth)
        ax.plot(esgm[i],rg_mean[i],'.',color=rgb,markersize=5*wth)
    ax.plot([],[],label='Capture Radius\n= %s $\\sigma_{ij}$'%csgm)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('Exclusion Radius ($\\sigma_{ij}$)',fontsize=size)
    ax.set_ylabel('Radius of Gyration ($\\mathrm{\\AA}$)',fontsize=size)
    ax.legend(loc='center',fontsize=size,handlelength=0.0)

    fig.subplots_adjust(right=figrig)
    plt.tight_layout(rect=[0,0,figrig,1])
    top=ax[-1].get_position().y1
    bot=ax[-1].get_position().y0
    cbar_ax=fig.add_axes([figrig+cbarwth,bot,cbarwth,top-bot])
    sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm,cax=cbar_ax)
    cbar.set_label('Exclusion Radius ($\\sigma_{ij}$)',fontsize=size)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(0.1))
    # cbar.set_ticks(esgm)
    cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
    cbar.outline.set_linewidth(wth)

    plt.savefig(f'{figname}_{csgm}sgm.png',format='png')
    plt.savefig(f'{figname}_{csgm}sgm.pdf',format='pdf')
    plt.show()