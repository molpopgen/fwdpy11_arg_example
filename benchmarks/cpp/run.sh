#!/bin/bash

# 72 hour run limit
SECONDS_TO_KILL=`echo "72*60*60"|bc -l`

for N in 1000 10000 100000
do
    for size in 1000 10000 100000
    do
        #We will arbitrarily GC every 0.1N generations
        GC=`echo "0.1*$N"|bc -l| sed 's/\.0//'`
        GC=1000
        TIME_MEM_FILE=time_arg.N$N"."size$size".out"
        DETAILED_TIME_FILE=detailed_time_arg.N$N"."size$size".out"
        /usr/bin/time -f "%e %M" -o $TIME_MEM_FILE parallel --timeout $SECONDS_TO_KILL PYTHONPATH=../.. python ../../benchmarking.py --popsize $N --theta $size --rho $size --nsam 100 --pdel 0.01 --gc $GC --seed ::: 42 > $DETAILED_TIME_FILE &
    done
done
