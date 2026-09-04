import numpy as np

nhi=1
dt=1

P=np.array([[[1]],[[1]],[[0]],[[1]],[[1]],[[0]],[[0]],[[1]],[[1]],[[0]]])
nfr,nch=np.shape(P)[:2]

Ppad=np.pad(P,((1,1),(0,0),(0,0)),'constant',constant_values=-1)
dP=Ppad[1:]-Ppad[:-1]
Ps=np.where(dP[:-1]!=0,1,0)
Pe=np.where(dP[1:]!=0,1,0)

bin=np.linspace(-1,3,nhi+1)
Len0=[]
Len1=[]
Len0log=[]
Len1log=[]
for i in range(nch):
    for j in range(i+1):
        seglen=[]
        segval=[]
        idxs=np.where(Ps[:,i,j]==1)[0]
        idxe=np.where(Pe[:,i,j]==1)[0]
        for k,l in zip(idxs,idxe):
            seglen+=[l-k+1]
            segval+=[P[k,i,j]]
        seglen=np.array(seglen)
        segval=np.array(segval)
        len0=seglen[segval==0]
        len1=seglen[segval==1]
        Len0+=(len0*dt).tolist()
        Len1+=(len1*dt).tolist()
        Len0log+=np.log10(np.array([x for x in len0 if x!=0])*dt).tolist()
        Len1log+=np.log10(np.array([x for x in len1 if x!=0])*dt).tolist()
print(Len0)
print(Len1)