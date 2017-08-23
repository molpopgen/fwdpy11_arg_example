import fwdpy11_arg_example.brute_force as bf
import sys

N = int(sys.argv[1])
rho = float(sys.argv[2])
seed = int(sys.argv[3])
timings = bf.evolve_track_wrapper(popsize=N, rho=rho, seed = seed)
print(timings)
