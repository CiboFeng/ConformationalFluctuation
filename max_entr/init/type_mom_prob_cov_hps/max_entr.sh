#!/bin/bash

itrs=0
itre=1
sysname=$(echo $(find . -type f -name "exp_cont_prob*") | cut -c 17- | cut -d '.' -f 1)
nmom=$(sed -n '2p' "$(find . -type f -name "exp_cont_prob*" -print -quit)" | awk '{print NF}')
njob=10

touch ./init.sh
echo -e "#!/bin/bash
#SBATCH --job-name init
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_init.out
#SBATCH --error ./job_out/%j_init.err

echo Begin" > ./init.sh
jobid=$(sbatch --parsable ./init.sh)
echo init.sh: $jobid

for i in $(seq $itrs 1 $itre)
do

	if [ $i == 0 ] 
	then
		mkdir ./sim_0
		mkdir ./coeff_0
		touch ./coeff_0/max_entr_init.py.sh
		echo -e "#!/bin/bash
#SBATCH --job-name init
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_init.out
#SBATCH --error ./job_out/%j_init.err

python ./max_entr_init.py -s ./seq_$sysname.pyw -c ./coeff_0/coeff_$sysname.pyw -d ./sim_0/sim_$sysname -nm $nmom -np $njob" > ./coeff_0/max_entr_init.py.sh
    	jobid=$(sbatch --dependency=afterok:$jobid --parsable ./coeff_0/max_entr_init.py.sh)
		echo $i max_entr_init.py.sh: $jobid
	fi

	mkdir ./sim_$(($i+1))
	mkdir ./pc_$i
	mkdir ./coeff_$(($i+1))

	touch ./coeff_$i/max_entr_file.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_file.out
#SBATCH --error ./job_out/%j_file.err

python ./max_entr_file.py -s ./seq_$sysname.pyw -co ./coeff_$i/coeff_$sysname.pyw -cf ./sim_$i/sim_$sysname -nf ./sim_$(($i+1))/sim_$sysname -np $njob" > ./coeff_$i/max_entr_file.py.sh
    jobid=$(sbatch --dependency=afterok:$jobid --parsable ./coeff_$i/max_entr_file.py.sh)
	echo $i max_entr_file.py.sh: $jobid

	touch ./sim_$i/max_entr_sim.sh
	echo -e "#!/bin/bash
#SBATCH --job-name sim
#SBATCH --ntasks 4
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%A_%a_sim.out
#SBATCH --error ./job_out/%A_%a_sim.err
#SBATCH --array 0-$(($njob-1))

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH
export PATH=/hpc2hdd/home/chu-amat/cbfengphy/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH
export PATH=/hpc2hdd/home/chuo657/cbfengphy/opt/lammps/install/lammps_09.10.2020_custom/bin:\$PATH

i=\$SLURM_ARRAY_TASK_ID
mpirun -n 4 lmp -in ./sim_$i/sim_${sysname}_\$i.cfg" > ./sim_$i/max_entr_sim.sh
	jobid=$(sbatch --dependency=afterok:$jobid --parsable ./sim_$i/max_entr_sim.sh)
	echo $i max_entr_sim.sh: $jobid
		
	touch ./pc_$i/max_entr_pc.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name pc
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%A_%a_pc.out
#SBATCH --error ./job_out/%A_%a_pc.err
#SBATCH --array 0-$(($njob-1))

i=\$SLURM_ARRAY_TASK_ID
python ./max_entr_pc.py -s ./seq_$sysname.pyw -t ./sim_$i/sim_${sysname}_\$i.trj -p ./pc_$i/pc_${sysname}_\$i.npz -n $nmom" > ./pc_$i/max_entr_pc.py.sh
	jobid=$(sbatch --dependency=afterany:$jobid --parsable ./pc_$i/max_entr_pc.py.sh)
	echo $i max_entr_pc.py.sh: $jobid

	touch ./coeff_$i/max_entr_coeff.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name coeff
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_coeff.out
#SBATCH --error ./job_out/%j_coeff.err

python ./max_entr_coeff.py -i ./coeff_$i/coeff_$sysname.pyw -d ./pc_$i -se ./seq_$sysname.pyw -p ./exp_cont_prob_$sysname.pyw -c ./chk_$sysname.pyw -o ./coeff_$(($i+1))/coeff_$sysname.pyw -sy $sysname" > ./coeff_$i/max_entr_coeff.py.sh
	jobid=$(sbatch --dependency=afterany:$jobid --parsable ./coeff_$i/max_entr_coeff.py.sh)
    echo $i max_entr_coeff.py.sh: $jobid

done
