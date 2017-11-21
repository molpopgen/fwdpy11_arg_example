#!/bin/bash

LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=medium.prof PYTHONPATH=../.. python ../../../benchmarking.py --popsize 50000 --theta 10000 --rho 10000 --nsam 100 --pdel 0.01 --gc 1000 --outfile1 medium.gz --seed 42
