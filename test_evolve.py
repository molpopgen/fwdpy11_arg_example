import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys
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

count = 0
for idx, pos in enumerate(evolver.sites.position):
    if(idx > 0 and pos == evolver.sites.position[idx-1]):
      #print(idx-1,pos,evolver.mutations[idx-1].node,evolver.unpack_index(evolver.mutations[idx-1].metadata), evolver.pop.mutations[evolver.unpack_index(evolver.mutations[idx-1].metadata)])
      #print(idx,pos,evolver.mutations[idx].node,evolver.unpack_index(evolver.mutations[idx].metadata), evolver.pop.mutations[evolver.unpack_index(evolver.mutations[idx].metadata)])
      count += 1

print(count)

count = 0
for mut in evolver.mutations:
    if(evolver.sites.position[mut.site] != evolver.pop.mutations[evolver.unpack_index(mut.metadata)].pos):
       #print(evolver.sites.position[mut.site], mut.site, mut.node, evolver.pop.mutations[evolver.unpack_index(mut.metadata)])
       count += 1
print(count)

msp_rng = msprime.RandomGenerator(seed+2)
neutral_sites = msprime.SiteTable()
neutral_mutations = msprime.MutationTable()
mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
mutgen.generate(evolver.nodes, evolver.edges, neutral_sites, neutral_mutations)


trees_selected = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=evolver.sites, mutations=evolver.mutations)
trees_neutral = msprime.load_tables(nodes=evolver.nodes, edges=evolver.edges, sites=neutral_sites, mutations=neutral_mutations)

home = str(Path.home())
for counter, tree in enumerate(trees_selected.trees()):
   if(tree.num_mutations > 0):
      tree.draw(path=home+"/tree_selected.svg", width=1500, height=1000, format="svg")
      break

for counter2, tree in enumerate(trees_neutral.trees()):
    if(counter2 == counter):
       tree.draw(path=home+"/tree_neutral.svg", width=1500, height=1000, format="svg")

