"""Import Modules"""
import numpy as np

""""""
dir='../../cont/cont_slab'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
outname='top_show'
nch=100
fr=0

""""""
P=np.zeros((nmd,nT,nch,nch))
for i in range(nmd):
    for j in range(nT):
        P[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')['Pchn'][fr]

""""""
np.save(f'{outname}.npy',P)