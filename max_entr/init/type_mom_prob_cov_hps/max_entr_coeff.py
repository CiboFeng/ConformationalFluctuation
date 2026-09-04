"""Import Modules"""
import argparse
import os
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/opt/mypylib')
sys.path.append('/hpc2hdd/home/chuo657/cbfengphy/opt/mypylib')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputting coefficient file')
parser.add_argument('-d',type=str,help='The directory of .npz files')
parser.add_argument('-se',type=str,help='The sequence file')
parser.add_argument('-p',type=str,help='The experimental contact probability file')
parser.add_argument('-c',type=str,help='The check file')
parser.add_argument('-o',type=str,help='The outputting coefficient file')
parser.add_argument('-sy',type=str,help='The system name')
args=parser.parse_args()
icoefffile=args.i
dir=args.d
seqfile=args.se
pexpfile=args.p
chkfile=args.c
ocoefffile=args.o
sysname=args.sy
ipathls=icoefffile.split('/')[:-1]
ipath=''
for i in range(len(ipathls)):
    ipath+=ipathls[i]+'/'
psimfile=ipath+f'sim_cont_prob_{sysname}.pyw'
ntp=20
ord=int(dir.split('_')[-1])
ord_1=20
scl_0=0.01

"""Read Sequence, Experimental Contact Probabilities and the Coefficients of the Last Step from Files"""
seq=pyw(seqfile,'Sequence','int')[0]
p_exp_all=pyw(pexpfile,'Probability')
coeff_i_all=pyw(icoefffile,'Coefficient')
idx=np.tril_indices(ntp)
N=len(seq)
nmo=len(p_exp_all)

P_exp=np.zeros((N,N,nmo))
for i in range(N):
    for j in range(N):
        if abs(i-j)>=4:
            P_exp[i,j]=p_exp_all[:,i*N+j]
p_exp=np.zeros((ntp,ntp,nmo))
for i in range(N):
    for j in range(i+1):
        p_exp[seq[i],seq[j],:]+=P_exp[i,j,:]
        if seq[i]!=seq[j]:
            p_exp[seq[j],seq[i],:]+=P_exp[i,j,:]
p_exp=p_exp[idx].T.reshape(-1)
Coeff_i=np.zeros((ntp,ntp,nmo))
for i in range(ntp):
    for j in range(ntp):
        Coeff_i[i,j]=coeff_i_all[:,i*ntp+j]
coeff_i=Coeff_i[idx]

"""Read Contact Probability Moment Matrics from .npz Files"""
files=[f for f in os.listdir(dir) if f.endswith('.npz')]
print('The Number of Trajectories:',len(files))
nfrs=np.zeros(len(files))
p_sim=np.zeros(round(ntp*(ntp+1)/2)*nmo)
C=np.zeros((round(ntp*(ntp+1)/2)*nmo,round(ntp*(ntp+1)/2)*nmo))
P_sim=np.zeros((N,N,nmo))
for i in range(len(files)):
    data=np.load(dir+f'/{files[i]}')
    nfrs[i]=data['nfr']
    p_sim+=nfrs[i]*data['p']
    C+=nfrs[i]*data['C']
    P_sim+=nfrs[i]*data['P']
    os.remove(dir+f'/{files[i]}')
print('The Number of Frames:',round(np.sum(nfrs)))
p_sim/=np.sum(nfrs)
C/=np.sum(nfrs)
C=C-np.dot(p_sim.reshape(-1,1),p_sim.reshape(1,-1))
P_sim/=np.sum(nfrs)

"""Update Coefficients"""
Cinv=np.linalg.pinv(C)
coeff_dlt_tmp=np.matmul(Cinv,(p_sim-p_exp).reshape(-1,1)).reshape(-1)
coeff_dlt=np.zeros((round(ntp*(ntp+1)/2),nmo))
for i in range(nmo):
    coeff_dlt[:,i]=coeff_dlt_tmp[i*round(ntp*(ntp+1)/2):(i+1)*round(ntp*(ntp+1)/2)]
scl=min(scl_0**((ord_1-ord)/ord_1),1)
print('Scaling Factor:',scl)
coeff_dlt*=scl
coeff_o=coeff_i+coeff_dlt

Coeff_o=np.zeros((ntp,ntp,nmo))
for i in range(len(idx[0])):
    Coeff_o[idx[0][i],idx[1][i]]=coeff_o[i]
    Coeff_o[idx[1][i],idx[0][i]]=coeff_o[i]

"""Compare Simulation Results with Experiment and Record"""
rel_err=[]
rel_err_log=[]
for i in range(nmo):
    rel_err.append([])
    rel_err_log.append([])
    for j in range(N):
        for k in range(N):
            if P_exp[j,k,i]+P_sim[j,k,i]!=0:
                rel_err[-1]+=[abs((P_sim[j,k,i]-P_exp[j,k,i])/(P_sim[j,k,i]+P_exp[j,k,i]))]
            if P_exp[j,k,i]!=0 and P_sim[j,k,i]!=0:
                rel_err_log[-1]+=[abs(np.log(P_sim[j,k,i]/P_exp[j,k,i])/np.log(P_sim[j,k,i]*P_exp[j,k,i]))]
rel_err=[np.mean(rel_err[i]) for i in range(nmo)]
rel_err_log=[np.mean(rel_err_log[i]) for i in range(nmo)]

p_exp_all=P_exp.reshape(-1,nmo)
p_sim_all=P_sim.reshape(-1,nmo)
r2=[r2_score(p_exp_all[:,i],p_sim_all[:,i]) for i in range(nmo)]
p_exp_all=[]
p_sim_all=[]
for i in range(nmo):
    p_exp_all.append([])
    p_sim_all.append([])
    for j in range(N):
        for k in range(N):
            if P_exp[j,k,i]!=0 and P_sim[j,k,i]!=0:
                p_exp_all[-1]+=[P_exp[j,k,i]]
                p_sim_all[-1]+=[P_sim[j,k,i]]
r2_log=[r2_score(np.log(p_exp_all[i]),np.log(p_sim_all[i])) for i in range(nmo)]

"""Output Data"""
if not os.path.exists(chkfile):
    chk=open(chkfile,'w')
    chk.write(f'The Number of Iteration Steps, the Relative Error of Contact Probabilities of 1-{nmo} Moments, Similar for the Logarithmic Relative Error, the Coefficient of Determination, and the Logarithmic Coefficient of Determination:\n')
    chk.close()
chk=open(chkfile,'a')
chk.write(f'{ord}\t')
for i in range(nmo):
    chk.write(f'{tot_len(rel_err[i],8)}\t')
for i in range(nmo):
    chk.write(f'{tot_len(rel_err_log[i],8)}\t')
for i in range(nmo):
    chk.write(f'{tot_len(r2[i],8)}\t')
for i in range(nmo):
    chk.write(f'{tot_len(r2_log[i],8)}\t')
chk.write('\n')
chk.close()

pyw=open(ocoefffile,'w')
pyw.write(f'The 1-{nmo}th Coefficient of Data-driven Term:\n')
for i in range(ntp):
    for j in range(ntp):
        for k in range(nmo):
            pyw.write(f'{Coeff_o[i,j,k]}\t')
        pyw.write('\n')
pyw.close()

pyw=open(psimfile,'w')
pyw.write(f'The 1-{nmo}th Order Contact Probability Moment:\n')
for i in range(N):
    for j in range(N):
        for k in range(nmo):
            pyw.write(f'{P_sim[i,j,k]}\t')
        pyw.write('\n')
pyw.close()

"""Plot"""
for i in range(nmo):
    plt.figure(figsize=figaspect(1/3))
    ax=plt.subplot(121)
    norm=colors.Normalize()
    plt.imshow(Coeff_i[:,:,i],norm=norm)
    plt.gca().invert_yaxis()
    plt.xlabel('Index of Particle Types')
    plt.ylabel('Index of Particle Types')
    cbar=plt.colorbar()
    cbar.set_label(f'The {i+1}-th Coefficient')
    
    ax=plt.subplot(122)
    norm=colors.SymLogNorm(0.001)
    plt.imshow(Coeff_i[:,:,i],norm=norm)
    plt.gca().invert_yaxis()
    plt.xlabel('Index of Particle Types')
    plt.ylabel('Index of Particle Types')
    cbar=plt.colorbar()
    cbar.set_label(f'The {i+1}-th Coefficient')
    plt.tight_layout()
    plt.savefig(ipath+f'coeff_{i+1}_{sysname}.png',format='png')
    plt.savefig(ipath+f'coeff_{i+1}_{sysname}.pdf',format='pdf')

Compare=np.zeros((N,N,nmo))
for i in range(N):
    for j in range(N):
        if i>j:
            for k in range(nmo):
                Compare[i,j,k]=P_exp[i,j,k]
        else:
            for k in range(nmo):
                Compare[i,j,k]=P_sim[i,j,k]
Compare[Compare==0]=np.min(Compare[Compare>0])/2

for i in range(nmo):
    plt.figure(figsize=figaspect(1/3))
    ax=plt.subplot(121)
    norm=colors.Normalize()
    plt.imshow(Compare[:,:,i],norm=norm)
    plt.gca().invert_yaxis()
    plt.title('Upper Left: Exp, Lower Right: Sim')
    plt.xlabel('Index of Particles')
    plt.ylabel('Index of Particles')
    cbar=plt.colorbar()
    cbar.set_label(f'The {i+1}th Order\nContact Probability Moment')
    
    ax=plt.subplot(122)
    norm=colors.LogNorm()
    plt.imshow(Compare[:,:,i],norm=norm)
    plt.gca().invert_yaxis()
    plt.title('Upper Left: Exp, Lower Right: Sim')
    plt.xlabel('Index of Particles')
    plt.ylabel('Index of Particles')
    cbar=plt.colorbar()
    cbar.set_label(f'The {i+1}th Order\nContact Probability Moment')
    plt.tight_layout()
    plt.savefig(ipath+f'cont_prob_{i+1}_{sysname}.png',format='png')
    plt.savefig(ipath+f'cont_prob_{i+1}_{sysname}.pdf',format='pdf')

for i in range(nmo):
    plt.figure(figsize=figaspect(1/3))
    ax=plt.subplot(121)
    plt.plot(p_exp_all[i],p_sim_all[i],'.')
    plt.plot([min(p_exp_all[i]),max(p_exp_all[i])],[min(p_exp_all[i]),max(p_exp_all[i])])
    plt.xlabel(f'The {i+1}th Order Contact\nProbability Moments of Experiment')
    plt.ylabel(f'The {i+1}th Order Contact\nProbability Moments of Simulation')
    
    ax=plt.subplot(122)
    plt.plot(p_exp_all[i],p_sim_all[i],'.')
    plt.plot([min(p_exp_all[i]),max(p_exp_all[i])],[min(p_exp_all[i]),max(p_exp_all[i])])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(f'The {i+1}th Order Contact\nProbability Moments of Experiment')
    plt.ylabel(f'The {i+1}th Order Contact\nProbability Moments of Simulation')
    plt.tight_layout()
    plt.savefig(ipath+f'cont_prob_cor_{i+1}_{sysname}.png',format='png')
    plt.savefig(ipath+f'cont_prob_cor_{i+1}_{sysname}.pdf',format='pdf')
