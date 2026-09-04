#!/bin/bash

for i in 1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm
do
	for j in 350.0 400.0 450.0
	do
		mkdir fus_${i}_${j}
		rm -r fus_${i}_${j}/job_out
		cd fus_${i}_${j}
		for k in $(seq 0 1 199)
		do 
			cp ../../fus_$i/sim_20/sim_fus_$k.cfg fus_$k.cfg
			sed -i 's#./sim_20/sim_fus_'$k'#fus_'$k'#g' fus_$k.cfg
                	sed -i 's#./sim_20/sim_fus#fus#g' fus_$k.cfg
                	sed -i 's#300.0#'$j'#g' fus_$k.cfg
			sed -i 's#custom 10000 fus_'$k'.trj id type x y z#xtc 10000 fus_'$k'.xtc#g' fus_$k.cfg
                        sed -i '/^dump_modify/d' fus_$k.cfg ### Can't use # as a separator.
			del_lines=1
        		tot_lines=$(wc -l < fus_$k.cfg)
        		let kp_lines=$tot_lines-$del_lines+1
        		sed -i "${kp_lines},\$d" fus_$k.cfg ### Can't use ''.
			cp ../../fus_$i/sim_20/sim_fus_$k.dat fus_$k.dat
		done
		cp ../../fus_$i/sim_20/sim_fus.par fus.par
		touch fus.sh
		echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --ntasks 4
#SBATCH --exclude cpu3-4,cpu3-5,cpu3-9,cpu3-12
#SBATCH --partition a128m512u
#SBATCH --output ./job_out/%A_%a.out
#SBATCH --error ./job_out/%A_%a.err
#SBATCH --array 0-199

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

i=\$SLURM_ARRAY_TASK_ID
mpirun -n 4 lmp -in fus_\$i.cfg" > fus.sh
		sbatch fus.sh
		cd ..
	done
done

