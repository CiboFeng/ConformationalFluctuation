#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_file.out
#SBATCH --error ./job_out/%j_file.err

python ./max_entr_file.py -s ./seq_fus.pyw -co ./coeff_14/coeff_fus.pyw -cf ./sim_14/sim_fus -nf ./sim_15/sim_fus -np 200
