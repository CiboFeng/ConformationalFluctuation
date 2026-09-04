#!/bin/bash

#dirs=(1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm)
#Ts=(350.0 400.0 450.0)
dirs=(1.5sgm_0.82sgm 1.55sgm_1.0sgm)
Ts=(500.0)

for i in $(seq 0 1 2)
do
	for j in $(seq 0 1 2)
	do
		sed -i "s|^coefffile=.*|coefffile='../../max_entr_all_atom/fus_"${dirs[$i]}"/coeff_20/coeff_fus.pyw'|"  lmp_file_decom.py
        	sed -i "s|^T=.*|T="${Ts[$j]}"|"  lmp_file_decom.py
		python lmp_file_decom.py
        	mkdir -p ../slab_${dirs[$i]}_${Ts[$j]}/rst
        	rm -r ../slab_${dirs[$i]}_${Ts[$j]}/job_out
        	mv fus.cfg fus.par fus.dat ../slab_${dirs[$i]}_${Ts[$j]}
        	cd ../slab_${dirs[$i]}_${Ts[$j]}
        	touch fus.sh
        	echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 128
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --exclude cpu3-3,cpu3-4,cpu3-5,cpu3-9,cpu3-12
#SBATCH --output ./job_out/%j.out
#SBATCH --error ./job_out/%j.err

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

mpirun -n 128 lmp -in fus.cfg" > fus.sh
		sbatch fus.sh
        	cd ../init
	done
done
