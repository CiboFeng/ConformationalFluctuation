#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_file.out
#SBATCH --error ./job_out/%j_file.err

python ./max_entr_file.py -s tdp43 -w . -co ./coeff_18/coeff_tdp43.pyw -cf ./tdp43.cfg -p ./lmp_18/tdp43.par -thermo 0 -Th 1000 -Tl 300 -na 500000 -nr 500000 -dump 2000 -ne 2000000
