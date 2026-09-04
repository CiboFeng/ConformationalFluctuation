"""Import Modules"""
import numpy as np
import random
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/opt/mypylib')
sys.path.append('/hpc2hdd/home/chuo657/cbfengphy/opt/mypylib')
from tot_len import tot_len

"""Set Arguments"""
sysname='fus'
T0=150
Tc=300
T1=800
NT=32
T='${T}'
nthermo=0
ndump=10000
nstep=1000000000
ntemper=1000
Kb=10
Rb=3.82
eps=0.2
lbox=1000
res1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
res3=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
lmd=[0.602942,0.558824,0.588236,0.294119,0.64706,0.558824,0.0,0.57353,0.764707,0.705883,0.720589,0.382354,0.676471,0.82353,0.758824,0.588236,0.588236,1.0,0.897059,0.664707]
culcut=np.zeros((20,20))
for i in [1,3,6,11]:
    for j in [1,3,6,11]:
        culcut[i,j]=35.0
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY ' \
    'GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN ' \
    'PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
N=len(seq)
seq=[res1.index(seq[i]) for i in range(N)]
ntp=len(res1)

"""Write to .cfg File"""
cfg=open(f'{sysname}.cfg','w')
cfg.write('units\t\treal\n')
cfg.write('boundary\tm m m\n')
cfg.write('timestep\t10.0\n')
cfg.write('dielectric\t\t80.0\n')
cfg.write('neighbor\t\t3.5 multi\n')
cfg.write('neigh_modify\tevery 10 delay 0\n')
cfg.write('special_bonds\tlj/coul 0.0 1.0 1.0\n\n')
cfg.write('atom_style\tfull\n')
cfg.write('bond_style\tharmonic\n')
cfg.write(f'pair_style\tlj/lambda 0.1 0.0 35.0\n\n')
cfg.write(f'read_data\t{sysname}.dat\n')
cfg.write(f'include\t\t{sysname}.par\n\n')
cfg.write(f'thermo\t\t{nthermo}\n')
cfg.write('thermo_style\tcustom step temp ke pe ebond epair\n\n')
cfg.write('minimize\t1.0e-6\t1.0e-8\t1000\t100000\n')
cfg.write('velocity\tall create 300.0 5725 rot yes dist gaussian\n')
cfg.write('fix\t\tbalance all balance 1000 1.1 shift xyz 10 1.01\n\n')
cfg.write('fix\t\tnve all nve\n\n')
cfg.write(f'variable\tI world {" ".join([str(int(x)) for x in np.linspace(1,NT,NT)])}\n')
Ts0=T0*(T1/T0)**((np.linspace(1,NT,NT)-1)/(NT-1))
N1=sum(1 for x in Ts0 if x<Tc)
fidx=np.log(Tc/T0)/np.log(T1/T0)*(NT-1)+1
N1+=int(np.sign(np.abs(N1-1-fidx+1)-np.abs(N1-fidx+1)))
Ts1=T0*(Tc/T0)**((np.linspace(1,N1,N1)-1)/(N1-1))
N2=NT-N1+1
Ts2=Tc*(T1/Tc)**((np.linspace(1,N2,N2)-1)/(N2-1))
Ts=np.concatenate((Ts1,Ts2[1:]))
Ts=[f'{float(x):.3f}' for x in Ts]
cfg.write(f'variable\tT world {" ".join(Ts)}\n')
cfg.write(f'fix\t\tlangevin all langevin {T} {T} 1000.0 5725 zero yes\n')
cfg.write(f'dump\t\tdump all custom {ndump} {sysname}_$I.trj id type x y z\n')
cfg.write('dump_modify\tdump sort id\n')
cfg.write(f'temper\t{nstep} {ntemper} {T} langevin 5725 5725\n\n')
cfg.close()

"""Write to .par File"""
par=open(f'{sysname}.par','w')
for i in range(ntp):
    par.write(f'mass\t{i+1}\t{tot_len(m[i],8)}\n')
par.write('\n')
par.write(f'bond_coeff\t1\t{tot_len(Kb,8)}\t{tot_len(Rb,8)}\n\n')
for i in range(ntp):
    for j in range(i,ntp):
        par.write(f'pair_coeff\t{i+1}\t{j+1}\t{tot_len(eps,8)}\t{tot_len((sgm[i]+sgm[j])/2,8)}\t'
                  f'{tot_len((lmd[i]+lmd[j])/2-0.08,8)}\t{tot_len(4*(sgm[i]+sgm[j])/2,8)}\t{tot_len(culcut[i,j],8)}\n')

"""Define the Structure of Protein"""
x=[0]
y=[0]
z=[0]
for j in range(1,N):
    judge=True
    while judge:
        theta=random.uniform(0,np.pi)
        phi=random.uniform(0,2*np.pi)
        x1=x[j-1]+Rb*np.sin(theta)*np.cos(phi)
        y1=y[j-1]+Rb*np.sin(theta)*np.sin(phi)
        z1=z[j-1]+Rb*np.cos(theta)
        dist=[]
        for k in range(len(x)):
            dist+=[((x1-x[k])**2+(y1-y[k])**2+(z1-z[k])**2)**0.5]
        dmin=min(dist)
        if dmin>Rb:
            judge=False
    x+=[x1]
    y+=[y1]
    z+=[z1]

"""Write to .dat File"""
dat=open(f'{sysname}.dat','w')
dat.write('LAMMPS Data\n\n')
dat.write(f'{N}\tatoms\n')
dat.write(f'{N-1}\tbonds\n\n')
dat.write(f'{ntp}\tatom types\n')
dat.write(f'1\tbond types\n\n')
dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\txlo\txhi\n')
dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\tylo\tyhi\n')
dat.write(f'{tot_len(-lbox/2,8)}\t{tot_len(lbox/2,8)}\tzlo\tzhi\n\n')
dat.write('Atoms\n\n')
for j in range(N):
    dat.write(f'{j+1}\t1\t{seq[j]+1}\t{tot_len(chrg[seq[j]],8)}\t{tot_len(x[j],8)}\t{tot_len(y[j],8)}\t{tot_len(z[j],8)}\n')
dat.write('\n')
dat.write('Bonds\n\n')
for j in range(N-1):
    dat.write(f'{j+1}\t1\t{j+1}\t{j+2}\n')
dat.write('\n')
dat.close()

