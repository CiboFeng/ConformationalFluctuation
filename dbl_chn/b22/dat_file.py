"""Import Modules"""
import argparse
import numpy as np
import os
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The .log files')
args=parser.parse_args()
logfile=args.f
mdl='_'.join(logfile.split('/')[-2].split('_')[-2:])
d0=logfile.split('/')[-1].split('_')[1]
outname=f'mass_cent_dist/{mdl}_{d0}'
S=0.8
nrep=3

"""Read Log Data"""
d=[]
for i in range(nrep):
    d_tmp=[]
    f=open(logfile.replace('_.',f'_{i}.'),'r')
    ls=f.readlines()
    for j in range(len(ls)):
        l_s=ls[j].strip('\n').split()
        if len(l_s)==7:
            if l_s[0]=='Step' and l_s[6]=='v_dist':
                k=j+1
                while True:
                    try:
                        l_s=ls[k].strip('\n').split()
                        d_tmp+=[float(l_s[6])]
                        k+=1
                    except Exception:
                        break
                break
    d+=d_tmp[-round(len(d_tmp)*S):]
nfr=len(d)
d=np.array(d)
z=np.abs((d-np.median(d))/(np.percentile(d,75)-np.percentile(d,25)))
d=d[z<10]

"""Output"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except Exception:
            pass
f=open(f'{outname}.dat','w')
for i in range(nfr):
    f.write(f'{i} {d[i]}\n')
f.close()
