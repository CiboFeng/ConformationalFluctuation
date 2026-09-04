#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_file.out
#SBATCH --error ./job_out/%j_file.err

python ./max_entr_file.py -s ./seq_fus.pyw -co ./coeff_30/coeff_fus.pyw -cf ./sim_30/sim_fus -nf ./sim_31/sim_fus -np 200
