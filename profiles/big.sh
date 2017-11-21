#!/bin/bash

LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=big.prof PYTHONPATH=../.. python ../../../benchmarking.py --popsize 50000 --theta 100000 --rho 100000 --nsam 100 --pdel 0.01 --gc 1000 --outfile1 big.gz --seed 42
