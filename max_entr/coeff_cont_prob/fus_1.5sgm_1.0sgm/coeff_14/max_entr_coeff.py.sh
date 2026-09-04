#!/bin/bash
#SBATCH --job-name coeff
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_coeff.out
#SBATCH --error ./job_out/%j_coeff.err

python ./max_entr_coeff.py -i ./coeff_14/coeff_fus.pyw -d ./pc_14 -se ./seq_fus.pyw -p ./exp_cont_prob_fus.pyw -c ./chk_fus.pyw -o ./coeff_15/coeff_fus.pyw -sy fus
