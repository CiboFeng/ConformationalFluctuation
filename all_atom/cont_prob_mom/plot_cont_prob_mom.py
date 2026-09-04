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
file='cont_prob_mom_1.5sgm.npy'
figname=file[:-4]
nmo=5

"""Read Data"""
P=np.load(file)
print(np.shape(P))
P=P.reshape(-1,nmo)

"""Fit Moments"""
def func(x,a):
    return x**a

Par=np.ones(nmo-1)
for i in range(1,nmo):
    par,cov=curve_fit(func,P[:,i-1],P[:,i])
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

"""Plot"""
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
lenbar=8

# fig,ax=plt.subplots(4,nmo,figsize=(4*nmo+3,15))
plt.figure(figsize=(4*(nmo-1)+3,10))
for i in range(1,nmo):
    ax=plt.subplot(3,nmo-1,i)
    ax.plot(P[:,i-1],P[:,i],'.',linewidth=wth)
    x=np.linspace(np.min(P[:,i-1]),np.max(P[:,i-1]),100)
    # ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
    a,b=par(i)
    ax.plot(x,a*x**b,linewidth=wth,label='$y=%.3fx^{%.3f}$'%(a,b))
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
    ax.set_xlabel('The %s-th Moment\nof Contact Probability'%(i),fontsize=size)
    ax.set_ylabel('The %s-th Moment\nof Contact Probability'%(i+1),fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend(loc='best',fontsize=size,handlelength=0.0)

    ax=plt.subplot(3,nmo-1,1*(nmo-1)+i)
    ax.plot(P[:,i-1],np.abs((a*P[:,i-1]**b-P[:,i])/P[:,i]),'.',linewidth=wth)
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
    ax.set_xlabel('The %s-th Moment\nof Contact Probability'%(i),fontsize=size)
    ax.set_ylabel('The Relative Error\nof %s-th Moment'%(i+1),fontsize=size)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax=plt.subplot(3,nmo-1,2*(nmo-1)+i)
    ax.plot(P[:,i-1],P[:,i],'.',linewidth=wth)
    P[:,i-1][P[:,i-1]==0]=np.nan
    x=np.logspace(np.log10(np.nanmin(P[:,i-1])),np.log10(np.nanmax(P[:,i-1])),100)
    # ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
    ax.plot(x,a*x**b,linewidth=wth,label='$y=%.3fx^{%.3f}$'%(a,b))
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
    ax.set_xlabel('The %s-th Moment\nof Contact Probability'%(i),fontsize=size)
    ax.set_ylabel('The %s-th Moment\nof Contact Probability'%(i+1),fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best',fontsize=size,handlelength=0.0)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

plt.figure(figsize=(5,10))
ax=plt.subplot(3,1,1)
ax.plot(P[:,0],P[:,1]-P[:,0]**2,'.',linewidth=wth)
x=np.linspace(np.min(P[:,0]),np.max(P[:,0]),100)
# ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
a,b=par(1)
ax.plot(x,a*x**b-x**2,linewidth=wth,label='$y=%.3fx^2$'%(a-1))
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
ax.set_xlabel('The Mean of Contact Probability',fontsize=size)
ax.set_ylabel('Variance',fontsize=size)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.legend(loc='best',fontsize=size,handlelength=0.0)

ax=plt.subplot(3,1,2)
ax.plot(P[:,0],np.abs((a*P[:,0]**b-P[:,1])/(P[:,1]-P[:,0]**2)),'.',linewidth=wth)
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
ax.set_xlabel('The Mean of Contact Probability',fontsize=size)
ax.set_ylabel('The Relative Error\nof Variance',fontsize=size)
# ax.set_xscale('log')
ax.set_yscale('log')

ax=plt.subplot(3,1,3)
ax.plot(P[:,0],P[:,1]-P[:,0]**2,'.',linewidth=wth)
P[:,0][P[:,0]==0]=np.nan
x=np.logspace(np.log10(np.nanmin(P[:,0])),np.log10(np.nanmax(P[:,0])),100)
# ax.plot(x/10**(2*i-1),x**Par[i-1]/10**(2*i+1),linewidth=wth,label='$y=x^{%.3f}$'%Par[i-1])
ax.plot(x,a*x**b-x**2,linewidth=wth,label='$y=%.3fx^2$'%(a-1))
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
ax.set_xlabel('The Mean of Contact Probability',fontsize=size)
ax.set_ylabel('Variance',fontsize=size)
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='best',fontsize=size,handlelength=0.0)

plt.tight_layout()
plt.savefig(f'{figname}_2.png',format='png')
plt.savefig(f'{figname}_2.pdf',format='pdf')
# plt.show()