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
file0='exp_cont_prob_fus.pyw'
file='coeff_20/sim_cont_prob_fus.pyw'
figname='cont_prob_cor'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
Mdls=['Qch','Mid','Flx']
nmd=len(mdls)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res_sort=['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res)
nat=len(seq)
res_=['A','N','D','Q','G','M','P','S','T','Y']
res_sort_=['D','S','T','N','Q','G','P','A','M','Y']
ntp_=len(res_)

"""Read Data"""
dseq=np.arange(nat)
xtp_=np.arange(ntp_)
num=np.zeros((ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        num[res.index(seq[i]),res.index(seq[j])]+=1
        if seq[i]!=seq[j]:
            num[res.index(seq[j]),res.index(seq[i])]+=1

P0=np.zeros((nmd,nat,nat))
P=np.zeros((nmd,nat,nat))
for i in range(nmd):
    P0_tmp=pyw(f'{dir}/fus_{mdls[i]}/{file0}','Probability')[0]
    P_tmp=pyw(f'{dir}/fus_{mdls[i]}/{file}','Probability')[0]
    for j in range(nat):
        for k in range(nat):
            P0[i,j,k]=P0_tmp[j*nat+k]
            P[i,j,k]=P_tmp[j*nat+k]

P_show=np.zeros((nmd,nat,nat))
for i in range(nat):
    for j in range(nat):
        if i>=j+4:
            P_show[:,i,j]=P0[:,i,j]
        if i<=j-4:
            P_show[:,i,j]=P[:,i,j]

Pl0=np.zeros((nmd,nat))
Pl=np.zeros((nmd,nat))
for i in range(nat):
    for j in range(nat-i):
        Pl0[:,i]+=P0[:,i+j,j]
        Pl[:,i]+=P[:,i+j,j]
Pl0/=np.arange(nat,0,-1)[np.newaxis]
Pl/=np.arange(nat,0,-1)[np.newaxis]

for i in range(nat):
    for j in range(nat):
        if abs(i-j)<4:
            P0[:,i,j]=0
            P[:,i,j]=0

p0=np.zeros((nmd,ntp,ntp))
p=np.zeros((nmd,ntp,ntp))
for i in range(nat):
    for j in range(i+1):
        p0[:,res.index(seq[i]),res.index(seq[j])]+=P0[:,i,j]
        p[:,res.index(seq[i]),res.index(seq[j])]+=P[:,i,j]
        if seq[i]!=seq[j]:
            p0[:,res.index(seq[j]),res.index(seq[i])]+=P0[:,i,j]
            p[:,res.index(seq[j]),res.index(seq[i])]+=P[:,i,j]
for i in range(ntp):
    for j in range(ntp):
        if num[i,j]!=0:
            p0[:,i,j]/=num[i,j]
            p[:,i,j]/=num[i,j]

idx_sort=[res.index(i) for i in res_sort]
p0=p0[:,idx_sort]
p0=p0[:,:,idx_sort]
p=p[:,idx_sort]
p=p[:,:,idx_sort]

p0_=np.zeros((nmd,ntp_,ntp_))
p_=np.zeros((nmd,ntp_,ntp_))
i_=0
for i in range(ntp):
    if res_sort[i] in res_:
        j_=0
        for j in range(ntp):
            if res_sort[j] in res_:
                p0_[:,i_,j_]=p0[:,i,j]
                p_[:,i_,j_]=p[:,i,j]
                j_+=1
        i_+=1

bin=50
p_show=np.zeros((nmd,ntp_*bin,ntp_*bin))
for i in range(ntp_*bin):
    for j in range(ntp_*bin):
        if i>j+1:
            p_show[:,i,j]=p0_[:,divmod(i,bin)[0],divmod(j,bin)[0]]
        if i<j+1:
            p_show[:,i,j]=p_[:,divmod(i,bin)[0],divmod(j,bin)[0]]
        if divmod(i,bin)[1]==0 or divmod(i,bin)[1]==bin-1 or divmod(j,bin)[1]==0 or divmod(j,bin)[1]==bin-1:
            p_show[:,i,j]=np.inf
        if abs(i-j)<=1:
            p_show[:,i,j]=np.inf

pl0=np.sum(p0_,axis=-1)
pl=np.sum(p_,axis=-1)

"""Plot"""
for i in range(nmd):
    fig=plt.figure(figsize=(25,15))
    gs=fig.add_gridspec(nrows=2,ncols=1,height_ratios=[1.3,1],hspace=0.3)
    Gs=[gs[0].subgridspec(nrows=1,ncols=4,width_ratios=[1,1,1,1],wspace=0.4),
        gs[1].subgridspec(nrows=1,ncols=4,width_ratios=[1,1,1,1],wspace=0.4)]
    axs=[[fig.add_subplot(Gs[i][j]) for j in range(Gs[i].ncols)] for i in range(len(Gs))]
    marglft=0.5
    margrig=0.5
    margbot=0.5
    margtop=0.5
    spcwth=0.3
    spchei=-0.1
    cbarwth=0.0
    cbarspcwth=0.0
    wth=3
    font={'family':'Arial','size':30}
    lenmaj=15
    lenmin=8
    lenbar=8
    xtick=0.1
    ytick=20
    ctick=1
    color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
    nodes=[0/3,1/3,2/3,3/3]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    ann=[['(A)','(C)','(E)','(G)'],['(B)','(D)','(F)','(H)']]

    ax=axs[0][0]
    ax.annotate(ann[0][0],xy=(0.0,1.3),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    norm=colors.Normalize(vmin=0,vmax=0.3)
    img=ax.imshow(P_show[i],cmap=cmap,norm=norm)
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
    ax.text(0.05,0.95,'AA',transform=ax.transAxes,va='top',ha='left',fontdict=font)
    ax.text(0.95,0.05,'CG',transform=ax.transAxes,va='bottom',ha='right',fontdict=font)
    cbar=plt.colorbar(img,orientation='horizontal',location='top')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.outline.set_linewidth(wth)
    cbar.set_label('$\\langle{}p_{ij}\\rangle$',fontdict=font)

    ax=axs[0][1]
    ax.annotate(ann[0][1],xy=(0.0,1.3),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    norm=colors.LogNorm(vmin=1e-4,vmax=1)
    img=ax.imshow(P_show[i],cmap=cmap,norm=norm)
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
    cbar.set_label('$\\langle{}p_{ij}\\rangle$',fontdict=font)

    ax=axs[0][2]
    ax.annotate(ann[0][2],xy=(0.0,1.3),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    norm=colors.Normalize(vmin=0,vmax=10)
    img=ax.imshow(p_show[i]*100,cmap=cmap,norm=norm)
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
    ax.set_xticks((xtp_+0.5)*bin,res_sort_,fontdict=font)
    ax.set_yticks((xtp_+0.5)*bin,res_sort_,fontdict=font)
    cbar=plt.colorbar(img,orientation='horizontal',location='top')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin/2,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    cbar.outline.set_linewidth(wth)
    cbar.ax.xaxis.set_major_locator(MultipleLocator(5))
    cbar.set_label('$\\langle{}p_{s_is_j}\\rangle\\times100$',fontdict=font)

    ax=axs[0][3]
    ax.annotate(ann[0][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.plot(dseq[4:],Pl0[i,4:],linewidth=wth,label='AA')
    ax.plot(dseq[4:],Pl[i,4:],linewidth=wth,label='CG')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$|i-j|$',fontdict=font)
    ax.set_ylabel('$\\langle{}p_{ij}\\rangle$',fontdict=font)
    ax.set_xscale('log')
    ax.set_yscale('log')
    leg=ax.legend(loc='best',handlelength=0.0,handletextpad=0.0,prop={'family':font['family'],'size':font['size']})
    for hnd,txt in zip(leg.legend_handles,leg.get_texts()):
        txt.set_color(hnd.get_color())

    ax=axs[1][0]
    ax.annotate(ann[1][0],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.plot(P0[i].reshape(-1),P[i].reshape(-1),'.',color='C2',linewidth=wth)
    x0=[0,max(np.max(P0[i]),np.max(P[i]))]
    ax.plot(x0,x0,color='k',linewidth=wth,label='$y=x$')
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\langle{}p_{ij}\\rangle$ of AA',fontdict=font)
    ax.set_ylabel('$\\langle{}p_{ij}\\rangle$ of CG',fontdict=font)
    ax.legend(loc='best',prop={'family':font['family'],'size':font['size']})

    ax=axs[1][1]
    ax.annotate(ann[1][1],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.plot(P0[i].reshape(-1),P[i].reshape(-1),'.',color='C2',linewidth=wth)
    x0=[min(np.min(P0[i][P0[i]!=0]),np.min(P[i][P[i]!=0])),max(np.max(P0[i]),np.max(P[i]))]
    ax.plot(x0,x0,color='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(0))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\langle{}p_{ij}\\rangle$ of AA',fontdict=font)
    ax.set_ylabel('$\\langle{}p_{ij}\\rangle$ of CG',fontdict=font)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax=axs[1][2]
    ax.annotate(ann[1][2],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    ax.scatter(p0[i].reshape(-1),p[i].reshape(-1),color='C2',s=20*wth)
    x0=[0,max(np.max(p0[i]),np.max(p[i]))]
    ax.plot(x0,x0,color='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$\\langle{}p_{s_is_j}\\rangle$ of AA',fontdict=font)
    ax.set_ylabel('$\\langle{}p_{s_is_j}\\rangle$ of CG',fontdict=font)

    ax=axs[1][3]
    ax.annotate(ann[1][3],xy=(0.0,1.05),xycoords='axes fraction',fontfamily=font['family'],fontsize=font['size'],fontweight='bold',va='bottom',ha='right')
    spc=0.4
    ax.bar(xtp_-spc/2,pl0[i],width=spc,color='C0',edgecolor='k',linewidth=0)
    ax.bar(xtp_+spc/2,pl[i],width=spc,color='C1',edgecolor='k',linewidth=0)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='x',which='both',direction='in',width=wth,length=0,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelfontfamily=font['family'],labelsize=font['size'])
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelfontfamily=font['family'],labelsize=font['size'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks(xtp_,res_sort_,fontdict=font)
    ax.set_xlabel('$s_i$',fontdict=font)
    ax.set_ylabel('$\\Sigma_{s_j}{}\\langle{}p_{s_is_j}\\rangle$',fontdict=font)

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

    plt.savefig(f'{figname}_{Mdls[i].lower()}.png',format='png',dpi=50)
    plt.savefig(f'{figname}_{Mdls[i].lower()}.pdf',format='pdf')
    # plt.show()