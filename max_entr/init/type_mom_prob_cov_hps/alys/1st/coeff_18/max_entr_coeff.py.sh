#!/bin/bash
#SBATCH --job-name coeff
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_coeff.out
#SBATCH --error ./job_out/%j_coeff.err

python ./max_entr_coeff.py -s tdp43 -w . -d ./pc_18 -i ./coeff_18/coeff_tdp43.pyw -o ./coeff_19/coeff_tdp43.pyw
