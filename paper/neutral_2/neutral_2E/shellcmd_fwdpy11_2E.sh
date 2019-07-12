#!/bin/bash

#$ -q krt2
#$ -pe openmp 105
cd $SGE_O_WORKDIR 

module load krthornt/anaconda
source activate fwdpy11_0_2_0

PYTHONPATH=/share/kevin2/dlawrie/fwdpy11_arg_example python ../../../test_evolve.py -o fwdpy11_sim2E.txt -g 2000 -1 flat 15000 -2 5000 1 1001 -nT 1000 -R 1000 -m 0 0 0.3333 3 3 0 -T 0 -B 0 -ns1 100 --anc_sam1 1000 5 --anc_sam2 1000 5 -S 14 -r 10000 > nohup_fwdpy11_2E.txt

