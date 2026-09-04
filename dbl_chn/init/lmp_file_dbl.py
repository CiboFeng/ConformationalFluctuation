"""Import Modules"""
import numpy as np
import random
import sys
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('D:\\Work\\Code\\Functions')
from tot_len import tot_len
from pyw import pyw
from xtc import xtc_rd
from perd_restr import str_restr

"""Set Arguments"""
coefffile='../../max_entr/fus_1.37sgm_0.6sgm/coeff_20/coeff_fus.pyw'
sysname='fus'
nch=2
nthermo=10000
T=300.0
ndump=10000
nstep=1000000000
nrestart=10000000
Kb=10
Rb=3.82
eps=0.2
lbox=500.0
mu=1.0
res1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res3=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
csgm=float(coefffile.split('/')[-3].split('_')[-2].replace('sgm',''))
esgm=float(coefffile.split('/')[-3].split('_')[-1].replace('sgm',''))
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY ' \
    'GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN ' \
    'PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
ntp=len(res1)
N=len(seq)
seq=[res1.index(seq[i]) for i in range(N)]
Kd=0.1
nd0=50
nrep=3
dlt_d0=3.4

"""Generate Equilibrium Distance"""
# d0=[0.0]
# for i in range(nd0-1):
#     step=3.4 if d0[i]<40.0 else 6.9
#     d0+=[d0[i]+step]
d0=np.round(np.arange(0.0,nd0*dlt_d0,dlt_d0),1)

"""Read Coefficient"""
coeff=pyw(coefffile,'Coefficient')[0]
Coeff=np.zeros((ntp,ntp))
for i in range(ntp):
    for j in range(ntp):
        Coeff[i,j]=coeff[i*ntp+j]

"""Write to .cfg File"""
for i in range(nd0):
    for j in range(nrep):
        cfg=open(f'{sysname}_{d0[i]}_{j}.cfg','w')
        cfg.write(f'log\t\t{sysname}_{d0[i]}_{j}.log\n')
        cfg.write('units\t\treal\n')
        cfg.write('boundary\tp p p\n')
        cfg.write('timestep\t10.0\n')
        cfg.write('dielectric\t\t80.0\n')
        cfg.write('neighbor\t\t3.5 multi\n')
        cfg.write('neigh_modify\tevery 10 delay 0\n')
        cfg.write('special_bonds\tlj/coul 0.0 1.0 1.0\n\n')
        cfg.write('atom_style\tfull\n')
        cfg.write('bond_style\tharmonic\n')
        cfg.write(f'pair_style\thybrid/overlay tanh/power 15.0 lj/lambda 0.1 0.0 35.0\n\n')
        cfg.write(f'read_data\t{sysname}_{d0[i]}_{j}.dat\n')
        cfg.write(f'include\t\t{sysname}.par\n\n')
        cfg.write('group\t\tmol1 molecule 1\n')
        cfg.write('group\t\tmol2 molecule 2\n')
        cfg.write(f'variable\tKd index {Kd}\n')
        cfg.write(f'variable\td0 index {d0[i]}\n\n')
        cfg.write(f'thermo\t\t{nthermo}\n')
        cfg.write('thermo_style\tcustom step temp ke pe ebond epair\n\n')
        cfg.write('minimize\t1.0e-6\t1.0e-8\t1000\t100000\n')
        cfg.write(f'velocity\tall create {T} 5725 rot yes dist gaussian\n')
        cfg.write('fix\t\tbalance all balance 1000 1.1 shift xyz 10 1.01\n\n')
        cfg.write('fix\t\tnve all nve\n')
        cfg.write(f'fix\t\tlangevin all langevin {T} {T} 1000000.0 5725\n\n')
        cfg.write('fix\t\tumbrella mol1 spring couple mol2 ${Kd} 0.0 0.0 0.0 ${d0}\n')
        cfg.write(f'variable\tdist equal v_d0+f_umbrella[4]/v_Kd\n')
        cfg.write(f'thermo\t\t{nthermo}\n')
        cfg.write('thermo_style\tcustom step temp ke pe ebond epair v_dist\n\n')
        cfg.write(f'dump\t\tdump all xtc {ndump} {sysname}_{d0[i]}_{j}.xtc\n')
        cfg.write(f'restart\t\t{nrestart} ./rst_{d0[i]}_{j}/{sysname}.rst\n')
        cfg.write(f'run\t\t{nstep}')
        cfg.close()

"""Write to .par File"""
par=open(f'{sysname}.par','w')
for i in range(ntp):
    par.write(f'mass\t{i+1}\t{tot_len(m[i],16)}\n')
par.write('\n')
par.write(f'bond_coeff\t1\t{tot_len(Kb,16)}\t{tot_len(Rb,16)}\n\n')
for i in range(ntp):
    for j in range(i+1):
        par.write(f'pair_coeff\t{j+1}\t{i+1}\ttanh/power\t{tot_len(Coeff[i,j],16)}\t{tot_len(mu,16)}\t'
                  f'{tot_len(csgm*(sgm[i]+sgm[j])/2,16)}\t{tot_len(1,16)}\t{tot_len(3*csgm*(sgm[i]+sgm[j])/2,16)}\n') ### I J..., I<J
        par.write(f'pair_coeff\t{j+1}\t{i+1}\tlj/lambda\t{tot_len(eps,16)}\t{tot_len(esgm*(sgm[i]+sgm[j])/2,16)}\t'
                  f'{tot_len(0.0,16)}\t{tot_len(4*esgm*(sgm[i]+sgm[j])/2,16)}\t{tot_len(np.where(chrg[i]*chrg[j]==0,0.0,35.0),16)}\n')
### If the pair style is used multiple times in the pair_style command, then an additional numeric
### argument must also be specified which is a number from 1 to M where M is the number of times
### the sub-style was listed in the pair style command. The extra number indicates which instance
### of the sub-style these coefficients apply to.
par.close()

for i in range(nd0):
    for j in range(nrep):
        """Generate an Initial Conformation"""
        x=[0]
        y=[0]
        z=[0]
        for k in range(1,N):
            judge=True
            while judge:
                theta=random.uniform(0,np.pi)
                phi=random.uniform(0,2*np.pi)
                x_tmp=x[k-1]+Rb*np.sin(theta)*np.cos(phi)
                y_tmp=y[k-1]+Rb*np.sin(theta)*np.sin(phi)
                z_tmp=z[k-1]+Rb*np.cos(theta)
                dist=[]
                for l in range(len(x)):
                    dist+=[((x_tmp-x[l])**2+(y_tmp-y[l])**2+(z_tmp-z[l])**2)**0.5]
                dmin=min(dist)
                if dmin>Rb:
                    judge=False
            x+=[x_tmp]
            y+=[y_tmp]
            z+=[z_tmp]

        """Copy the Proteins to Multiple Positions"""
        if i==0:
            xc=[-d0[1]/2,d0[1]/2]
        else:
            xc=[-d0[i]/2,d0[i]/2]
        yc=[0,0]
        zc=[0,0]

        """Write to .dat File"""
        dat=open(f'{sysname}_{d0[i]}_{j}.dat','w')
        dat.write('LAMMPS Data\n\n')
        dat.write(f'{N*nch}\tatoms\n')
        dat.write(f'{(N-1)*nch}\tbonds\n\n')
        dat.write(f'{ntp}\tatom types\n')
        dat.write(f'1\tbond types\n\n')
        dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\txlo\txhi\n')
        dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\tylo\tyhi\n')
        dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\tzlo\tzhi\n\n')
        dat.write('Atoms\n\n')
        for k in range(nch):
            for l in range(N):
                dat.write(f'{k*N+l+1}\t{k+1}\t{seq[l]+1}\t{tot_len(chrg[seq[l]],8)}\t{tot_len(x[l]+xc[k],8)}\t{tot_len(y[l]+yc[k],8)}\t{tot_len(z[l]+zc[k],8)}\n')
        dat.write('\n')
        dat.write('Bonds\n\n')
        for k in range(nch):
            for l in range(N-1):
                dat.write(f'{k*(N-1)+l+1}\t1\t{k*N+l+1}\t{k*N+l+2}\n')
        dat.write('\n')
        dat.close()