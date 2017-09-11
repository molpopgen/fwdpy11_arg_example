import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys

N = int(sys.argv[1])
rho = float(sys.argv[2])
theta = float(sys.argv[3])
gc_interval = int(sys.argv[4])
seed = int(sys.argv[5])
simplifier,atracker,tsim = ea.evolve_track_wrapper(popsize=N, rho=rho, seed = seed, gc_interval=gc_interval, mu = 0.0)

print(tsim,simplifier.times)
np.random.seed(seed)

# Get a sample of size n = 10 
msprime.simplify_tables(np.random.choice(2*N, 10, replace = False).tolist(), nodes = simplifier.nodes, edgesets = simplifier.edgesets)
msp_rng = msprime.RandomGenerator(seed)
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
mutgen.generate(simplifier.nodes, simplifier.edgesets, sites, mutations)
print(sites.num_rows)
