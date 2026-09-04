"""Import Modules"""
import os
import numpy as np

"""Set Arguments"""
datdir='/hpc2hdd/home/cfeng593/proj/alys/rg_var/dbl_chn/b22/mass_cent_dist'
outdir='cfg'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Kd=0.1
nd0=22

"""Generate Equilibrium Distance"""
d0=[0.0]
for i in range(nd0-1):
    step=3.4 if d0[i]<40.0 else 6.9
    d0+=[d0[i]+step]
    
""""""
if not os.path.exists(f'{outdir}'):
    try:
        os.makedirs(outdir)
    except Exception:
        pass
for i in range(nmd):
    cfg=open(f'{outdir}/{mdls[i]}.cfg','w')
    for j in range(nd0):
        cfg.write(f'{datdir}/{mdls[i]}_{j}.dat {d0[j]} {Kd}\n')
    cfg.close()