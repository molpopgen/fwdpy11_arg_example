#!/bin/bash

# 72 hour run limit
SECONDS_TO_KILL=259200

export PYTHONPATH=../.. 
for N in 1000 10000 100000
do
    for size in 1000 10000 100000
    do
        #We will arbitrarily GC every 0.1N generations
        GC=1000
        TIME_MEM_FILE=time_with_neutral.N$N"."size$size".out"
        if [ ! -e $TIME_MEM_FILE ]
        then
            echo "SECONDS_TO_KILL = $SECONDS_TO_KILL"
            /usr/bin/time -f "%e %M" -o $TIME_MEM_FILE timeout $SECONDS_TO_KILL python ../../benchmarking.py --neutral_mutations --popsize $N --theta $size --rho $size --nsam 100 --pdel 0.01 --gc $GC --seed 42 &
        fi
    done
done
