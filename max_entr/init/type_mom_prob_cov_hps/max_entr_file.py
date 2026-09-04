"""Import Modules"""
import argparse
import numpy as np
import random
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/opt/mypylib')
sys.path.append('/hpc2hdd/home/chuo657/cbfengphy/opt/mypylib')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-s',type=str,help='The sequence file')
parser.add_argument('-co',type=str,help='The inputting coefficient file')
parser.add_argument('-cf',type=str,help='The directory and name of current files')
parser.add_argument('-nf',type=str,help='The directory and name of next files')
parser.add_argument('-np',type=int,help='The number of parallel trajectories')
args=parser.parse_args()
seqfile=args.s
coefffile=args.co
curpre=args.cf
nxtpre=args.nf
n=args.np
ntp=20
nthermo=0
Thi=1000.0
Tlo=300.0
nstep_annl=1000000
nstep_relx=1000000
ndump=5000
nstep_equl=5000000
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
Kb=10.0
R0=3.82
mu=1.0
csgm=1.5
esgm=1.0
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
eps=0.2
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]
dspace=10.0

"""Read the Sequence File and Coefficient File"""
seq=pyw(seqfile,'Sequence','int')[0]
N=len(seq)
coeff=pyw(coefffile,'Coefficient')
nmo=len(coeff)
Coeff=np.zeros((ntp,ntp,nmo))
for i in range(ntp):
    for j in range(ntp):
        Coeff[i,j]=coeff[:,i*ntp+j]

"""Write to .cfg Files Cyclically"""
for i in range(0,n):
    cfg=open(f'{curpre}_{i}.cfg','w')
    cfg.write(f'log\t\t{curpre}_{i}.log\n')
    cfg.write('units\t\treal\n')
    cfg.write('boundary\tm m m\n')
    cfg.write('timestep\t10.0\n')
    cfg.write('dielectric\t\t80.0\n')
    cfg.write('neighbor\t\t3.5 multi\n')
    cfg.write('neigh_modify\tevery 10 delay 0\n')
    cfg.write('special_bonds\tlj/coul 0.0 1.0 1.0\n\n')
    cfg.write('atom_style\tfull\n')
    cfg.write('bond_style\tharmonic\n')
    ps=''
    for j in range(nmo):
        ps+='tanh/power 15.0 '
    cfg.write(f'pair_style\thybrid/overlay {ps}lj/lambda 0.1 0.0 35.0\n\n')
    cfg.write(f'read_data\t{curpre}_{i}.dat\n')
    cfg.write(f'include\t\t{curpre}.par\n\n')
    cfg.write(f'thermo\t\t{nthermo}\n')
    cfg.write('thermo_style\tcustom step temp ke pe ebond epair\n\n')
    cfg.write('minimize\t1.0e-6\t1.0e-8\t1000\t100000\n')
    cfg.write(f'velocity\tall create {Thi} 5725 rot yes dist gaussian\n')
    cfg.write('fix\t\tbalance all balance 1000 1.1 shift xyz 10 1.01\n\n')
    cfg.write('fix\t\tnve all nve\n\n')
    cfg.write(f'fix\t\tlangevin1 all langevin {Thi} {Tlo} 1000.0 5725 zero yes\n')
    cfg.write(f'run\t\t{nstep_annl}\n')
    cfg.write('unfix\t\tlangevin1\n\n')
    cfg.write(f'fix\t\tlangevin2 all langevin {Tlo} {Tlo} 1000.0 5725 zero yes\n')
    cfg.write(f'run\t\t{nstep_relx}\n')
    cfg.write('unfix\t\tlangevin2\n\n')
    cfg.write(f'fix\t\tlangevin3 all langevin {Tlo} {Tlo} 1000.0 5725 zero yes\n')
    cfg.write(f'dump\t\tdump all custom {ndump} {curpre}_{i}.trj id type x y z\n')
    cfg.write('dump_modify\tdump sort id\n')
    cfg.write(f'run\t\t{nstep_equl}\n\n')
    cfg.write(f'write_data\t{nxtpre}_{i}.dat')
    cfg.close()

"""Write to .par File"""
par=open(f'{curpre}.par','w')
for i in range(ntp):
    par.write(f'mass\t{i+1}\t{tot_len(m[i],16)}\n')
par.write('\n')
par.write(f'bond_coeff\t1\t{tot_len(Kb,16)}\t{tot_len(R0,16)}\n\n')
for i in range(ntp):
    for j in range(i+1):
        if nmo==1:
            par.write(f'pair_coeff\t{j+1}\t{i+1}\ttanh/power\t{tot_len(Coeff[i,j,0],16)}\t{tot_len(mu,16)}\t'
                      f'{tot_len(csgm*(sgm[i]+sgm[j])/2,16)}\t{tot_len(1,16)}\t{tot_len(3*csgm*(sgm[i]+sgm[j])/2,16)}\n') ### I J..., I<J
        else:
            for k in range(nmo):
                par.write(f'pair_coeff\t{j+1}\t{i+1}\ttanh/power\t{k+1}\t{tot_len(Coeff[i,j,k],16)}\t{tot_len(mu,16)}\t'
                          f'{tot_len(csgm*(sgm[i]+sgm[j])/2,16)}\t{tot_len(k+1,16)}\t{tot_len(3*csgm*(sgm[i]+sgm[j])/2,16)}\n') ### I J..., I<J
        par.write(f'pair_coeff\t{j+1}\t{i+1}\tlj/lambda\t{tot_len(eps,16)}\t{tot_len(esgm*(sgm[i]+sgm[j])/2,16)}\t'
                  f'{tot_len(0.0,16)}\t{tot_len(4*esgm*(sgm[i]+sgm[j])/2,16)}\t{tot_len(np.where(chrg[i]*chrg[j]==0,0.0,35.0),16)}\n')
### If the pair style is used multiple times in the pair_style command, then an additional numeric
### argument must also be specified which is a number from 1 to M where M is the number of times
### the sub-style was listed in the pair style command. The extra number indicates which instance
### of the sub-style these coefficients apply to.
par.close()

"""Write to .dat Files Cyclically"""
for i in range(0,n):
    #"""Generate an Initial Conformation"""
    x=[0]
    y=[0]
    z=[0]
    for j in range(1,N):
        judge=True
        while judge:
            theta=random.uniform(0,np.pi)
            phi=random.uniform(0,2*np.pi)
            x1=x[j-1]+R0*np.sin(theta)*np.cos(phi)
            y1=y[j-1]+R0*np.sin(theta)*np.sin(phi)
            z1=z[j-1]+R0*np.cos(theta)
            dist=[]
            for k in range(len(x)):
                dist+=[((x1-x[k])**2+(y1-y[k])**2+(z1-z[k])**2)**0.5]
            dmin=min(dist)
            if dmin>R0:
                judge=False
        x+=[x1]
        y+=[y1]
        z+=[z1]

    #"""Calculate the Spacial Range of the Box"""
    xlo=min(x)-dspace*(max(x)-min(x))
    xhi=max(x)+dspace*(max(x)-min(x))
    ylo=min(x)-dspace*(max(y)-min(y))
    yhi=max(x)+dspace*(max(y)-min(y))
    zlo=min(x)-dspace*(max(z)-min(z))
    zhi=max(x)+dspace*(max(z)-min(z))

    #"""Write to .dat File"""
    dat=open(f'{nxtpre}_{i}.dat','w')
    dat.write('LAMMPS Data\n\n')
    dat.write(f'{N}\tatoms\n')
    dat.write(f'{N-1}\tbonds\n\n')
    dat.write(f'{ntp}\tatom types\n')
    dat.write(f'1\tbond types\n\n')
    dat.write(f'{tot_len(xlo,16)}\t{tot_len(xhi,16)}\txlo\txhi\n')
    dat.write(f'{tot_len(ylo,16)}\t{tot_len(yhi,16)}\tylo\tyhi\n')
    dat.write(f'{tot_len(zlo,16)}\t{tot_len(zhi,16)}\tzlo\tzhi\n\n')
    dat.write('Atoms\n\n')
    for j in range(N):
        dat.write(f'{j+1}\t1\t{seq[j]+1}\t{tot_len(chrg[seq[j]],16)}\t{tot_len(x[j],16)}\t{tot_len(y[j],16)}\t{tot_len(z[j],16)}\n')
    dat.write('\n')
    dat.write('Bonds\n\n')
    for j in range(N-1):
        dat.write(f'{j+1}\t1\t{j+1}\t{j+2}\n')
    dat.write('\n')
    dat.close()