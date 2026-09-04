#!/bin/bash

dirs=(1.37sgm_0.6sgm 1.5sgm_0.82sgm 1.55sgm_1.0sgm)
nodes=(cpu3-8 cpu3-13 cpu3-16)
nodes_excl=cpu3-1,cpu3-3,cpu3-4,cpu3-5,cpu3-6,cpu3-7,cpu3-9,cpu3-11,cpu3-12,cpu3-18,cpu3-19,cpu3-20

for i in $(seq 0 1 2)
do
        sed -i "s|^coefffile=.*|coefffile='../../max_entr/fus_"${dirs[$i]}"/coeff_20/coeff_fus.pyw'|"  lmp_file_sph.py
        python lmp_file_sph.py
        mkdir -p ../sph_${dirs[$i]}/rst
        rm -r ../sph_${dirs[$i]}/job_out
        mv fus.cfg fus.par fus.dat ../sph_${dirs[$i]}
        cd ../sph_${dirs[$i]}
        touch fus.sh
        echo -e "#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 128
#SBATCH --exclude ${nodes_excl}
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
