#!/bin/env sh

for gc in 10 19 73 100 200 10000000000000
do 
    PYTHONPATH=../../msprime python prototype_mutations_history_gc_debug.py --gc $gc 
done
