"""Import Modules"""
import os
import numpy as np

"""Set Arguments"""
datdir='/hpc2hdd/home/cfeng593/proj/alys/rg_var/dbl_chn/b22/mass_cent_dist_rep'
outdir='cfg_rep'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Kd=0.1
nd0=50
dlt_d0=3.4
nrep=3

"""Generate Equilibrium Distance"""
# d0=[0.0]
# for i in range(nd0-1):
#     step=3.4 if d0[i]<40.0 else 6.9
#     d0+=[d0[i]+step]
d0=np.round(np.arange(0.0,nd0*dlt_d0,dlt_d0),1)
    
""""""
if not os.path.exists(f'{outdir}'):
    try:
        os.makedirs(outdir)
    except Exception:
        pass
for i in range(nmd):
    for j in range(nrep):
        cfg=open(f'{outdir}/{mdls[i]}_{j}.cfg','w')
        for k in range(nd0):
            cfg.write(f'{datdir}/{mdls[i]}_{d0[k]}_{j}.dat {d0[k]} {Kd}\n')
        cfg.close()