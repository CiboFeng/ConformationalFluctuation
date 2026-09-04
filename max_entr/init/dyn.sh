#!/bin/bash

for i in 1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm
do
	mkdir fus_$i
	rm -r fus_$i/job_out
	cd fus_$i
	cp ../../fus_$i/sim_20/sim_fus_0.cfg fus.cfg
	cp ../../fus_$i/sim_20/sim_fus_0.dat fus.dat
	cp ../../fus_$i/sim_20/sim_fus.par fus.par
	touch fus.sh
	echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --ntasks 4
#SBATCH --exclude cpu3-4,cpu3-5,cpu3-9,cpu3-12
#SBATCH --partition a128m512u
#SBATCH --output ./job_out/%j.out
#SBATCH --error ./job_out/%j.err

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

mpirun -n 4 lmp -in fus.cfg" > fus.sh
	sed -i 's#./sim_20/sim_fus_0#fus#g' fus.cfg
	sed -i 's#./sim_20/sim_fus#fus#g' fus.cfg
	del_lines=14
	tot_lines=$(wc -l < fus.cfg)
	let kp_lines=$tot_lines-$del_lines+1
	sed -i "${kp_lines},\$d" fus.cfg
	echo -e "fix\t\tlangevin all langevin 300.0 300.0 1000.0 5725 zero yes
dump\t\tdump all xtc 10000 fus.xtc
run\t\t100000000" >> fus.cfg
	sbatch fus.sh
	cd ..
done

