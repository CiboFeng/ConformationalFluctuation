import numpy as np

f=open('vmd_rg.txt')
lsf=f.readlines()
Rg=[]
for i in range(1,len(lsf)):
    Rg+=[float(lsf[i].strip('\n').split()[-1])]
Rg_sel=[]
Rg_sel.append([np.min(Rg),np.min(Rg)*2/3+np.max(Rg)*1/3,np.min(Rg)*1/3+np.max(Rg)*2/3,np.max(Rg)])
Rg_sel.append([Rg_sel[0][1],Rg_sel[0][1]*2/3+Rg_sel[0][2]*1/3,Rg_sel[0][1]*1/3+Rg_sel[0][2]*2/3,Rg_sel[0][2]])
Rg_sel.append([Rg_sel[1][1],Rg_sel[1][1]*2/3+Rg_sel[1][2]*1/3,Rg_sel[1][1]*1/3+Rg_sel[1][2]*2/3,Rg_sel[1][2]])
# print(Rg_sel)

ord_sel=[]
for i in range(len(Rg_sel)):
    ord_sel+=[[]]
    for j in range(len(Rg_sel[i])):
        ord_sel[i]+=[[]]
        d=np.abs(np.array(Rg)-Rg_sel[i][j]).tolist()
        ord_sel[i][j]=d.index(np.min(d))+1
print(ord_sel)
