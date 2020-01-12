#!/bin/bash

#$ -q krt2
#$ -pe openmp 105
cd $SGE_O_WORKDIR 

module load krthornt/anaconda
source activate fwdpy11_0_2_0

PYTHONPATH=/share/kevin2/dlawrie/fwdpy11_arg_example python ../../test_continuity.py 42 9 > nohup_continuity_9.txt
