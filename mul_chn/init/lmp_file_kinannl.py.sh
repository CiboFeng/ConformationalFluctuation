#!/bin/bash

dirs=(1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm)
nodes=(cpu3-8 cpu3-13 cpu3-16)
nrep=20

for i in $(seq 0 1 2)
do
	sed -i "s|^coefffile=.*|coefffile='../../max_entr_all_atom/fus_"${dirs[$i]}"/coeff_20/coeff_fus.pyw'|"  lmp_file_kinannl.py
	sed -i 's|^nrep=.*|nrep='$nrep'|'  lmp_file_kinannl.py
	python lmp_file_kinannl.py
	for j in $(seq 0 1 $(($nrep-1)))
	do
		mkdir -p ../kinannl_${dirs[$i]}/rst_$j
	done
	rm -r ../kinannl_${dirs[$i]}/job_out
	mv fus*.cfg fus.par fus*.dat ../kinannl_${dirs[$i]}
	cd ../kinannl_${dirs[$i]}
	touch fus.sh
	echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 64
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --exclude cpu1-1,cpu1-4,cpu1-7,cpu1-8,cpu1-15,cpu1-19,cpu1-21,cpu1-23,cpu1-24,cpu1-25,cpu1-27,cpu1-29,cpu1-30,cpu1-31,cpu1-33,cpu1-34,cpu1-35,cpu1-37,cpu1-38,cpu1-43,cpu1-44,cpu1-46,cpu1-48,cpu1-51,cpu1-55,cpu1-60,cpu1-61,cpu1-62,cpu1-63,cpu1-64,cpu1-65,cpu1-66,cpu1-67,cpu1-69,cpu1-70,cpu1-71,cpu1-72,cpu1-74,cpu1-78,cpu1-82,cpu1-83,cpu1-93,cpu1-95,cpu1-98,cpu1-100,cpu1-102,cpu1-103,cpu1-106,cpu1-107
#SBATCH --output ./job_out/%A_%a.out
#SBATCH --error ./job_out/%A_%a.err
#SBATCH --array 0-$(($nrep-1))

i=\$SLURM_ARRAY_TASK_ID

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

mpirun -n 64 lmp -in fus_\$i.cfg" > fus.sh
	sbatch fus.sh
	cd ../init
done
