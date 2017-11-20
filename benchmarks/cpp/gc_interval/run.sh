#!/bin/bash

# 72 hour run limit
SECONDS=`echo "72*60*60"|bc -l`

for N in 1000 10000
do
    for size in 1000 10000 50000
    do
        for GC in 10 100 1000 10000
        do
            TIME_MEM_FILE=time_arg.N$N"."size$size"."GC$GC".out"
            DETAILED_TIME_FILE=detailed_time_arg.N$N"."size$size"."GC$GC".out"
            /usr/bin/time -f "%e %M" -o $TIME_MEM_FILE parallel --timeout $SECONDS PYTHONPATH=../../.. python ../../../benchmarking.py --popsize $N --theta $size --rho $size --nsam 100 --pdel 0.01 --gc $GC --seed ::: 42 > $DETAILED_TIME_FILE
        done
    done
done
