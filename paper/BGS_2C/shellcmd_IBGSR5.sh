#!/bin/bash

#$ -q krt2
#$ -pe openmp 64-91

cd $SGE_O_WORKDIR 

module load krthornt/anaconda
source activate fwdpy11_0_2_0

PYTHONPATH=/share/kevin2/dlawrie/fwdpy11_arg_example python ../../test_evolve.py -o fwdpy11_sim2C_IBGS_R5.txt -g 2000 -1 flat 10000 -2 500 1 1001 -nT 1000 -R 5000 -m 0 0 0.05 3 3 0 -T 0 -B 200000 -niT -T 250 -s=-0.01 -ns1 100 --anc_sam1 1000 5 --anc_sam2 1000 5 -nc $CORES -S 22 -r 90 > nohup_IBGS_R5.txt
