import fwdpy11_arg_example.brute_force as bf
import numpy as np
import sys

N = int(sys.argv[1])
rho = float(sys.argv[2])
seed = int(sys.argv[3])
simplifier,atracker = bf.evolve_track_wrapper(popsize=N, rho=rho, seed = seed)

print(np.array(atracker.nodes,copy=False))
