mdl=0 ### Index of model.
tclfile=f'vmd_seg_tube_{mdl}.tcl'
nch=100
N=163

f=open(tclfile,'w')
for i in range(nch):
    f.write('mol color Name\n')
    f.write('mol representation Tube 0.300000 12.000000\n')
    f.write('mol selection all\n')
    f.write('mol material Opaque\n')
    f.write(f'mol addrep {mdl}\n')
    f.write(f'mol modselect {i+1} {mdl} index {i*N+2} to {(i+1)*N-3}\n')
    f.write(f'mol modcolor {i+1} {mdl} ColorID {divmod(i,33)[1]}\n')
f.write(f'mol delrep 0 {mdl}\n')
f.close()