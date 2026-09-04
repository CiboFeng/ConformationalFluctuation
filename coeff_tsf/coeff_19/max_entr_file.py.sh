#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition a128m512u
#SBATCH --time 00:10:00
#SBATCH --output ./job_out/%j_file.out
#SBATCH --error ./job_out/%j_file.err

python ./max_entr_file.py -s fus_wt -w . -co ./coeff_19/coeff_fus_wt.pyw -cf ./fus_wt.cfg -p ./lmp_19/fus_wt.par -thermo 0 -Th 1000 -Tl 300 -na 1000000 -nr 1000000 -dump 10000 -ne 10000000
