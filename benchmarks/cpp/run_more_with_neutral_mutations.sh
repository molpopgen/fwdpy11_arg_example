#!/bin/bash

# 72 hour run limit
SECONDS_TO_KILL=259200

export PYTHONPATH=../.. 

for N in 1000 10000 
do
    KILL=0
    for size in 2500 5000 7500 15000
    do
        #We will arbitrarily GC every 0.1N generations
        echo "PARAMS = " $N $size
        GC=1000
        TIME_MEM_FILE=time_with_neutral.N$N"."size$size".out"
        if [ $KILL == 0 ]
        then
            DETAILED_TIME_FILE=detailed_time_argi_with_neutral.N$N"."size$size".out.gz"
            /usr/bin/time -f "%e %M" -o $TIME_MEM_FILE parallel --timeout $SECONDS_TO_KILL PYTHONPATH=../.. python ../../benchmarking.py --popsize $N --neutral_mutations --theta $size --rho $size --nsam 100 --pdel 0.01 --gc $GC --outfile1 $DETAILED_TIME_FILE --seed ::: 42 
            # If we terminated the job,
            # then the next params will
            # also take too long, and 
            # so we can skip running them
            STATUS=$?
            if [ $STATUS -gt 0 ]
            then
                echo "STATUS = " $STATUS
                KILL=1
            fi
        fi
            #/usr/bin/time -f "%e %M" -o $TIME_MEM_FILE timeout $SECONDS_TO_KILL python ../../benchmarking.py --neutral_mutations --popsize $N --theta $size --rho $size --nsam 100 --pdel 0.01 --gc $GC --seed 42
    done
done
