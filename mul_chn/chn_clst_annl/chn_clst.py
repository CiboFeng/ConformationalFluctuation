"""Import Modules"""
import argparse
import numpy as np
import networkx as nx
import os
from scipy.spatial.distance import squareform,pdist
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The .npz file')
args=parser.parse_args()
npzfile=args.f
mdl='_'.join(npzfile.split('/')[-1].split('_')[:2])
rep=npzfile.split('/')[-1][:-4].split('_')[2]
outname=f'chn_clst/{mdl}_{rep}'

""""""
P=np.load(npzfile)['Pchn']
nfr,nch=np.shape(P)[:2]
P=np.heaviside(P-1.0,0.0)
clst=[]
for i in range(nfr):
    G=nx.from_numpy_array(P[i])
    clst.append([list(x) for x in list(nx.connected_components(G))])

# r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-8]}.dat')
# nfr,N,_=np.shape(r)
# nch=round(N/nat)
# Sgm=np.zeros(nat)
# for i in range(nat):
#     Sgm[i]=sgm[seq[i]]
# Sgm=np.tile(Sgm,nch)
# Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2
# clst=[]
# for i in range(nfr):
#     # P=0.5*(1-np.tanh(mu*(squareform(pdist(r[i],'euclidean'))-csgm*Sgm)))
#     P=1-np.heaviside(squareform(pdist(r[i],'euclidean'))-csgm*Sgm,0.0)
#     # p=np.zeros((nch,nch))
#     # for j in range(nch):
#     #     for k in range(nch):
#     #         p[j,k]=np.heaviside(np.sum(P[j*nat:(j+1)*nat,k*nat:(k+1)*nat])-1.0,1.0)
#     p=np.heaviside(np.sum(P.reshape(nch,nat,nch,nat),axis=(1,3))-1.0,1.0)
#     G=nx.from_numpy_array(p)
#     clst.append([list(x) for x in list(nx.connected_components(G))])

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Frame Index, Chain Index, Cluster Index:\n')
for i in range(nfr):
    for j in range(len(clst[i])):
        for k in range(len(clst[i][j])):
            pyw.write(f'{i} {clst[i][j][k]} {j}\n')

# import numpy as np
# import networkx as nx
# import matplotlib.pyplot as plt
# x=np.array([[1,0,0,1,0,0,0,0],
#             [0,1,1,0,1,0,1,0],
#             [0,1,1,0,0,0,1,0],
#             [1,0,0,1,0,1,0,0],
#             [0,1,0,0,1,0,1,0],
#             [0,0,0,1,0,1,0,0],
#             [0,1,1,0,1,0,1,0],
#             [0,0,0,0,0,0,0,1]])
# y=np.array([[1,1,0,0,0,0,0,0],
#             [1,1,1,0,0,0,0,0],
#             [0,1,1,0,0,0,0,0],
#             [0,0,0,1,1,1,1,0],
#             [0,0,0,1,1,0,1,0],
#             [0,0,0,1,0,1,1,0],
#             [0,0,0,1,1,1,1,0],
#             [0,0,0,0,0,0,0,1]])
# G=nx.from_numpy_array(x)
# cpnt=[list(x) for x in list(nx.connected_components(G))]
# # cpnt=sorted([sorted(list(x)) for x in list(nx.connected_components(G))],key=lambda x:x[0])
# print(cpnt)
# # x=np.sort(np.sort(x,axis=0),axis=1)
# # sort_idx=[0,3,5,1,2,4,6,7]
# # x=x[sort_idx][:,sort_idx]
# # print(x)
