"""Import Modules"""
import os
import numpy as np

"""Set Arguments"""
dir='mass_cent_dist'
outname='distr'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
nhi=20
nd0=50
dlt_d0=3.4

"""Generate Equilibrium Distance"""
# d0=[0.0]
# for i in range(nd0-1):
#     step=3.4 if d0[i]<40.0 else 6.9
#     d0+=[d0[i]+step]
d0=np.round(np.arange(0.0,nd0*dlt_d0,dlt_d0),1)

"""Read Data"""
D=np.zeros((nmd,nd0,nhi))
P=np.zeros((nmd,nd0,nhi))
for i in range(nmd):
    for j in range(nd0):
        d=[]
        f=open(f'{dir}/{mdls[i]}_{d0[j]}.dat','r')
        ls=f.readlines()
        for k in range(len(ls)):
            l_s=ls[k].strip('\n').split()
            d+=[float(l_s[-1])]
        hist,edge=np.histogram(d,bins=nhi)
        D[i,j]=(edge[:-1]+edge[1:])/2
        P[i,j]=hist/np.sum(hist*(np.max(D[i,j])-np.min(D[i,j]))/nhi)

"""Output"""
np.savez(f'{outname}.npz',D=D,P=P)