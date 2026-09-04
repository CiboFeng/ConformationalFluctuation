#!/bin/bash
pwd=$(pwd)
find $pwd -type f -name 'log.lammps' -exec rm -r {} +
find $pwd -type f -name 'bash.log' -exec rm -r {} +
find $pwd -type f -name 'init.sh' -exec rm -r {} +
find $pwd -type d -name 'job_out' -exec rm -r {} +
find $pwd -type d -name 'sim_*' -exec rm -r {} +
find $pwd -type d -name 'pc_*' -exec rm -r {} +
find $pwd -type d -name 'coeff_*' -exec rm -r {} +
find $pwd -type f -name 'chk_*' -exec rm -r {} +
find $pwd -type f -name '*cfg' -exec rm -r {} +
find $pwd -type f -name 'times.dat' -exec rm -r {} +
find $pwd -type f -name 'trn_pair_*' -exec rm -r {} +
