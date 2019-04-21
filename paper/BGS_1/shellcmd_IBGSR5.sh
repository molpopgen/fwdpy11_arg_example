#!/bin/bash

#$ -q krt2
#$ -pe openmp 64-105

cd $SGE_O_WORKDIR 

module load krthornt/anaconda
source activate fwdpy11_0_2_0

PYTHONPATH=/share/kevin2/dlawrie/fwdpy11_arg_example python ../../test_evolve.py -o fwdpy11_sim1_IBGS_R5.txt -1 flat 3000 -2 3000 1 501 -nT 1000 -R 1500 -m 0 0 0.05 3 3 0 -B 60000 -niT -T 75 -s=-0.01 -ns1 100 --anc_sam1 500 10 --anc_sam2 500 10 -nc $CORES -S 124 -r 104 > nohup_IBGS_R5.txt
