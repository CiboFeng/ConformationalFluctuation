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
file='dist.npz'
figname=file[:-4]
nat=163
nhist=50
nmo=5

"""Read Data"""
data=np.load(file)
D=data['D']
p_D=data['p_D']

"""Calculate Moments"""
dx=D[:,:,1]-D[:,:,0]
xm=np.zeros((nat,nat,nmo))
for i in range(nmo):
    xm[:,:,i]=np.sum(D**(i+1)*p_D,axis=-1)*dx
xm=xm.reshape(-1,nmo)

"""Fit Moments"""
def func(x,a):
    return x**a

Par=np.ones(nmo-1)
for i in range(1,nmo):
    par,cov=curve_fit(func,xm[:,i-1],xm[:,i])
    Par[i-1]=par[0]
    print(f'{i}:\t{cov}')

""""""
def par(n):
    import math

    # def factorial2(n):
    #     if n<0:
    #         return None
    #     return math.prod(range(n,0,-2))
    #
    # if divmod(n,2)[1]==0:
    #     m=int(n/2)
    #     k=(2*m+1)/(2*m)
    #     a=2**(m+2)*math.factorial(m+1)/(np.sqrt(2*np.pi)*factorial2(2*m+1)**k)
    #     b=k
    # else:
    #     m=int((n+1)/2)
    #     k=2*m/(2*m-1)
    #     a=np.sqrt(2*np.pi)**k*factorial2(2*m+1)/(2**((m+1)*k)*math.factorial(m)**k)
    #     b=k

    a=(np.pi/4)**(1/2/n)*math.gamma((n+4)/2)/math.gamma((n+3)/2)**((n+1)/n)
    b=(n+1)/n

    return a,b

dseq=np.arange(nat)

"""Plot"""
fig,ax=plt.subplots(1,1,figsize=(6,4))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
lenbar=8
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
cmap=plt.cm.rainbow
norm=colors.Normalize(vmin=min(dseq),vmax=max(dseq))

for i in range(int(nat/2)-1):
    for j in [nat-i-1]:
        normc=norm(dseq[abs(i-j)])
        rgb=cmap(normc)
        ax.plot(D[i,j],p_D[i,j],color=rgb,linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(-7,165)
ax.set_xlabel(f'Distance',fontsize=size)
ax.set_ylabel(f'Probability Density',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')

fig.subplots_adjust(right=0.77)
plt.tight_layout(rect=[0,0,0.77,1])
top=ax.get_position().y1
bot=ax.get_position().y0
cbar_ax=fig.add_axes([0.8,bot,0.03,top-bot])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance',fontsize=size)
# cbar.ax1.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax1.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1))) ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(labelsize=size,direction='in',width=wth,length=lenmin)
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

# fig,ax=plt.subplots(4,nmo,figsize=(4*nmo+3,15))
plt.figure(figsize=(4*(nmo-1)+3,10))
for i in range(1,nmo):
    ax=plt.subplot(3,nmo-1,i)
    ax.plot(xm[:,i-1]/10**(2*i-1),xm[:,i]/10**(2*i+1),'.',linewidth=wth)
    D=np.linspace(np.min(xm[:,i-1]),np.max(xm[:,i-1]),100)
    # ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
    a,b=par(i)
    ax.plot(D/10**(2*i-1),a*D**b/10**(2*i+1),linewidth=wth,label='$y=%.3fx^{%.3f}$'%(a,b))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(-7,165)
    ax.set_xlabel('The %s-th Moment\nof Distance $\\times{}10^%s$'%(i,2*i-1),fontsize=size)
    ax.set_ylabel('The %s-th Moment\nof Distance $\\times{}10^%s$'%(i+1,2*i+1),fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend(loc='best',fontsize=size,handlelength=0.0)

    ax=plt.subplot(3,nmo-1,1*(nmo-1)+i)
    ax.plot(xm[:,i-1]/10**(2*i-1),np.abs((a*xm[:,i-1]**b-xm[:,i])/xm[:,i]),'.',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(-7,165)
    ax.set_xlabel('The %s-th Moment\nof Distance $\\times{}10^%s$'%(i,2*i-1),fontsize=size)
    ax.set_ylabel('The Relative Error\nof %s-th Moment'%(i+1),fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax=plt.subplot(3,nmo-1,2*(nmo-1)+i)
    ax.plot(xm[:,i-1]/10**(2*i-1),xm[:,i]/10**(2*i+1),'.',linewidth=wth)
    xm[:,i-1][xm[:,i-1]==0]=np.nan
    D=np.logspace(np.log10(np.nanmin(xm[:,i-1])),np.log10(np.nanmax(xm[:,i-1])),100)
    # ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
    ax.plot(D/10**(2*i-1),a*D**b/10**(2*i+1),linewidth=wth,label='$y=%.3fx^{%.3f}$'%(a,b))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    # ax.set_xlim(-7,165)
    ax.set_xlabel('The %s-th Moment\nof Distance $\\times{}10^%s$'%(i,2*i-1),fontsize=size)
    ax.set_ylabel('The %s-th Moment\nof Distance $\\times{}10^%s$'%(i+1,2*i+1),fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best',fontsize=size,handlelength=0.0)

plt.tight_layout()
plt.savefig(f'{figname}_cor.png',format='png')
plt.savefig(f'{figname}_cor.pdf',format='pdf')
# plt.show()

plt.figure(figsize=(5,10))
ax=plt.subplot(3,1,1)
ax.plot(xm[:,0],xm[:,1]-xm[:,0]**2,'.',linewidth=wth)
D=np.linspace(np.nanmin(xm[:,0]),np.nanmax(xm[:,0]),100)
# ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
a,b=par(1)
ax.plot(D,a*D**b-D**2,linewidth=wth,label='$y=%.3fx^2$'%(a-1))
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(-7,165)
ax.set_xlabel('The Mean of Distance',fontsize=size)
ax.set_ylabel('Variance',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.legend(loc='best',fontsize=size,handlelength=0.0)

ax=plt.subplot(3,1,2)
ax.plot(xm[:,0],np.abs((a*xm[:,0]**b-xm[:,1])/(xm[:,1]-xm[:,0]**2)),'.',linewidth=wth)
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(-7,165)
ax.set_xlabel('The Mean of Distance',fontsize=size)
ax.set_ylabel('The Relative Error\nof Variance',fontsize=size)
# ax.set_xscale('log')
ax.set_yscale('log')

ax=plt.subplot(3,1,3)
ax.plot(xm[:,0],xm[:,1]-xm[:,0]**2,'.',linewidth=wth)
xm[:,0][xm[:,0]==0]=np.nan
D=np.logspace(np.log10(np.nanmin(xm[:,0])),np.log10(np.nanmax(xm[:,0])),100)
# ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
ax.plot(D,a*D**b-D**2,linewidth=wth,label='$y=%.3fx^2$'%(a-1))
ax.autoscale()
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ax.set_xlim(-7,165)
ax.set_xlabel('The Mean of Distance',fontsize=size)
ax.set_ylabel('Variance',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='best',fontsize=size,handlelength=0.0)

plt.tight_layout()
plt.savefig(f'{figname}_cor2.png',format='png')
plt.savefig(f'{figname}_cor2.pdf',format='pdf')
# plt.show()