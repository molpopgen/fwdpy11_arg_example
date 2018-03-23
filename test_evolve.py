import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys

N = int(sys.argv[1])
rho = float(sys.argv[2])
theta = float(sys.argv[3])
gc_interval = int(sys.argv[4])
seed = int(sys.argv[5])
simplifier,tsim = ea.evolve_track_wrapper(popsize=N, rho=rho, seed=seed, gc_interval=gc_interval, mu=theta/float(4*N))

print(tsim,simplifier.times)
np.random.seed(seed)

# Get a sample of size n = 10 
msprime.simplify_tables(np.random.choice(2*N, 10, replace = False).tolist(), nodes = simplifier.nodes, edges = simplifier.edges, sites = simplifier.sites, mutations = simplifier.mutations)
msp_rng = msprime.RandomGenerator(seed)
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
mutgen.generate(simplifier.nodes, simplifier.edges, sites, mutations)
print(sites.num_rows)
