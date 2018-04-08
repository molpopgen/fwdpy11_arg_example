import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys
import pickle
from pathlib import Path

N = int(sys.argv[1])
rho = float(sys.argv[2])
theta = float(sys.argv[3])
gc_interval = int(sys.argv[4])
seed = int(sys.argv[5])

evolver = ea.evolve_track_wrapper(popsize=N, rho=rho, seed=seed, gc_interval=gc_interval, mu=theta/float(4*N))
print(evolver.times)

# Get a sample of size n = 10 
np.random.seed(seed+1)
curr_samples = np.random.choice(2*N, 10, replace = False).tolist()
samples = curr_samples+evolver.anc_samples
msprime.simplify_tables(samples, nodes = evolver.nodes, edges = evolver.edges, sites = evolver.sites, mutations = evolver.mutations)
print(evolver.sites.num_rows)

# count = 0
# for idx, pos in enumerate(evolver.sites.position):
#     if(idx > 0 and pos == evolver.sites.position[idx-1]):
#       for idx_mut in evolver.mutations.site:
#           if(idx_mut == idx or idx_mut == idx-1):
#              print(idx_mut,pos,evolver.mutations[idx_mut].node,pickle.loads(evolver.mutations[idx_mut].metadata))
#       count += 1
# 
# print(count)

msp_rng = msprime.RandomGenerator(seed+2)
neutral_sites = msprime.SiteTable()
neutral_mutations = msprime.MutationTable()
mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
mutgen.generate(evolver.nodes, evolver.edges, neutral_sites, neutral_mutations)


trees_selected = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=evolver.sites, mutations=evolver.mutations)
trees_neutral = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=neutral_sites, mutations=neutral_mutations)

home = str(Path.home())
trees_selected.first().draw(path=home+"/tree_selected.svg", width=1500, height=1000, format="svg")
trees_neutral.first().draw(path=home+"/tree_neutral.svg", width=1500, height=1000, format="svg")

