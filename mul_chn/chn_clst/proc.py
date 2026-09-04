"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import numpy as np

"""Set Arguments"""
dir='chn_clst'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
outname='proc'
nfr=10000
nrep=20
nch=50
dt=10e-6*10000 ### ns.
nsms=0.2 ### nsms order of magnitude.

"""Read Data and Calculate"""
t=np.arange(nfr)*dt
Q=np.zeros((nmd,nrep,nfr,nch),dtype=int)
nclst=np.zeros((nmd,nrep,nfr),dtype=int)
nisol=np.zeros((nmd,nrep,nfr),dtype=int)
nsq=np.zeros((nmd,nrep,nfr),dtype=int)
dnmin=np.zeros((nmd,nrep,nfr-1,round(nch/2)),dtype=int)
for i in range(nmd):
    for j in range(nrep):
        q=pyw(f'{dir}/{mdls[i]}_{j}.pyw','Cluster')
        nfr=round(np.shape(q)[1]/nch)
        for k in range(np.shape(q)[1]):
                Q[i,j,int(q[0,k]),int(q[1,k])]=int(q[2,k])
        for k in range(nfr):
            nclst[i,j,k]=len(list(set(Q[i,j,k])))
            num1=0
            num=0
            for l in range(np.max(Q[i,j,k])+1):
                num_=np.sum(Q[i,j,k]==l)
                if num_==1:
                    num1+=1
                num+=num_**2
            nisol[i,j,k]=num1
            nsq[i,j,k]=num
        for k in range(1,nfr):
            dnmin_=[]
            for l in range(np.max(Q[i,j,k])+1):
                Q_tmp=Q[i,j,k-1,Q[i,j,k]==l]
                num=[]
                for m in set(Q_tmp):
                    num+=[np.sum(Q_tmp==m)]
                if len(num)>1:
                    dnmin_+=[np.min(num)]
            dnmin[i,j,k-1,:len(dnmin_)]=dnmin_

treq=np.zeros((nmd,nrep,nch))
for i in range(nmd):
    for j in range(nrep):
        for k in range(nch):
            for l in range(nfr-1):
                if nclst[i,j,l]>k+1 and nclst[i,j,l+1]<=k+1:
                    treq[i,j,k]=((nclst[i,j,l+1]-k-1)*t[l]+(k+1-nclst[i,j,l])*t[l+1])/(nclst[i,j,l+1]-nclst[i,j,l])
                    break

dnsq=np.sqrt((nsq[:,:,1:]-nsq[:,:,:-1])/2)

nclst_mean=np.mean(nclst,axis=1)
nclst_std=np.std(nclst,axis=1)
treq_mean=np.mean(treq,axis=1)
treq_std=np.std(treq,axis=1)
nisol_mean=np.mean(nisol,axis=1)
nisol_std=np.std(nisol,axis=1)
dnsq_mean=np.mean(dnsq,axis=1)
dnsq_std=np.std(dnsq,axis=1)
dnmin_mean=np.mean(dnmin,axis=1)
dnmin_std=np.std(dnmin,axis=1)
# dnsq_mean=np.zeros((nmd,nfr-1))
# dnsq_std=np.zeros((nmd,nfr-1))
# dnmin_mean=np.zeros((nmd,nfr-1))
# dnmin_std=np.zeros((nmd,nfr-1))
# dnsq[dnsq==0.0]=np.nan
# dnmin=dnmin.astype(float)
# dnmin[dnmin==0.0]=np.nan
# for i in range(nfr-1):
#     hlf=nsms/2
#     idxl=round(np.floor(1/dt*10**(np.log10(i*dt)-hlf)))
#     idxh=round(np.ceil(1/dt*10**(np.log10(i*dt)+hlf)))
#     dnsq_mean[:,i]=np.nanmean(dnsq[:,:,idxl:min(idxh,nfr-1)],axis=(1,2))
#     dnsq_std[:,i]=np.nanstd(dnsq[:,:,idxl:min(idxh,nfr-1)],axis=(1,2))
#     dnmin_mean[:,i]=np.nanmean(dnmin[:,:,idxl:min(idxh,nfr-1)],axis=(1,2,3))
#     dnmin_std[:,i]=np.nanstd(dnmin[:,:,idxl:min(idxh,nfr-1)],axis=(1,2,3))

"""Output Data"""
np.savez(f'{outname}.npz',
         t=t,
         nclst=nclst,
         nclst_mean=nclst_mean,
         nclst_std=nclst_std,
         treq=treq,
         treq_mean=treq_mean,
         treq_std=treq_std,
         nisol=nisol,
         nisol_mean=nisol_mean,
         nisol_std=nisol_std,
         nsq=nsq,
         dnsq=dnsq,
         dnsq_mean=dnsq_mean,
         dnsq_std=dnsq_std,
         dnmin=dnmin,
         dnmin_mean=dnmin_mean,
         dnmin_std=dnmin_std)

# import numpy as np
# dt=0.1
# hlf=0.5
# for i in [1,2,3,30,300,3000,9999]:
#     idxl=round(1/dt*10**(np.log10(i*dt)-hlf))
#     idxh=round(1/dt*10**(np.log10(i*dt)+hlf))
#     print(idxl,min(idxh,9999))

# import numpy as np
# y=np.array([0,0,1,1,1,2,2,3,4])
# z=np.array([0,0,0,0,0,1,1,1,0])
# idx=[0,3,4,1,5,2,6,8,7]
# y=y[idx]
# z=z[idx]
# n=[]
# for i in range(np.max(z)+1):
#     yi=y[z==i]
#     n.append([])
#     for j in set(yi):
#         n[-1]+=[int(np.sum(yi==j))]
# print(n)