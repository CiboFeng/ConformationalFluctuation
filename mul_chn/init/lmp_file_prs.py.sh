#!/bin/bash

dirs=(1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm)
nodes=(cpu3-8 cpu3-13 cpu3-16)

for i in $(seq 0 1 2)
do
	sed -i "s|^coefffile=.*|coefffile='../../max_entr_all_atom/fus_"${dirs[$i]}"/coeff_20/coeff_fus.pyw'|"  lmp_file_prs.py
	python lmp_file_prs.py
	mkdir -p ../prs_${dirs[$i]}/rst
	rm -r ../prs_${dirs[$i]}/job_out
	mv fus.cfg fus.par fus.dat ../prs_${dirs[$i]}
	cd ../prs_${dirs[$i]}
	touch fus.sh
	echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 128
#SBATCH --nodelist ${nodes[$i]}
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --output ./job_out/%j.out
#SBATCH --error ./job_out/%j.err

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

mpirun -n 128 lmp -in fus.cfg" > fus.sh
	sbatch fus.sh
	cd ../init
done
