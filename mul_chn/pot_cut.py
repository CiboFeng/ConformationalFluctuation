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
figname='pot_cut'

""""""
def lj(x,eps=0.1,re=1.0):
    return 4*eps*((re/x)**12-(re/x)**6)
def soft(x,alf=-0.5,mu=7,rc=1.0):
    return alf/2*(1-np.tanh(mu*(x-rc)))
def debye(x,C=1.0,q1=1.0,q2=1.0,D=80.0,kpp=0.1):
    return C*q1*q2/D/x*np.exp(-kpp*x)

""""""
plt.figure(figsize=(8,5))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20

ax=plt.subplot(1,1,1)
x=np.arange(0.9,5.0,0.01)
ax.plot(x,np.abs(lj(x)),linewidth=wth,color='r',label='$U(r)=|4\\times{}0.2\\times{}0.5[(\\frac{1}{r})^{12}-(\\frac{1}{r})^6]|$')
ax.plot(x,np.abs(soft(x)),linewidth=wth,color='g',label='$U(r)=|-0.5/2[1-\\tanh7(r-1)]|$')
ax.plot(x,debye(7*x),linewidth=wth,color='b',label='$U(r)=|\\frac{1}{80(7r)}\\mathrm{e}^{-0.1(7r)}|$')
ax.plot(x,1/80/(7*x),linewidth=wth,color=(0.5,0.5,0.5),label='$U(r)=|\\frac{1}{80r}|$')
plt.plot([4,4],[0,np.abs(lj(4))],linewidth=wth,color='r',linestyle=':')
plt.plot([0,4],[np.abs(lj(4)),np.abs(lj(4))],linewidth=wth,color='r',linestyle=':')
plt.plot([2,2],[0,np.abs(soft(2))],linewidth=wth,color='g',linestyle=':')
plt.plot([0,2],[np.abs(soft(2)),np.abs(soft(2))],linewidth=wth,color='g',linestyle=':')
plt.plot([3,3],[0,np.abs(debye(7*3))],linewidth=wth,color='b',linestyle=':')
plt.plot([0,3],[np.abs(debye(7*3)),np.abs(debye(7*3))],linewidth=wth,color='b',linestyle=':')
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
ax.set_xlim(0.7,5.2)
ax.set_ylim(1e-8,1e1)
ax.set_yscale('log')
ax.set_xlabel('$r$',fontsize=size)
ax.set_ylabel('$U(r)$ (kcal/mol)',fontsize=size)
ax.legend(loc='best',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()