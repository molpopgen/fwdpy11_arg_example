import fwdpy11_arg_example.evolve_arg as ea
import msprime
import numpy as np
import sys

N = int(sys.argv[1])
rho = float(sys.argv[2])
theta = float(sys.argv[3])
gc_interval = int(sys.argv[4])
seed = int(sys.argv[5])
simplifier,atracker = ea.evolve_track_wrapper(popsize=N, rho=rho, seed = seed, gc_interval=100, mu = 0.0)

np.random.seed(seed)

# Get a sample of size n = 10 
msprime.simplify_tables(np.random.choice(2*N, 10, replace = False).tolist(), nodes = simplifier.nodes, edgesets = simplifier.edgesets)
msp_rng = msprime.RandomGenerator(seed)
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
mutgen.generate(simplifier.nodes, simplifier.edgesets, sites, mutations)
print(sites.num_rows)
# This next set of steps is a little dumb:
# We are over-writing the data for our population,
# which is not desirable.
# nodes = simplifier.nodes
# edgesets = simplifier.edgesets
# x = msprime.load_tables(nodes=nodes, edgesets=edgesets)
# samples = np.random.choice(2*N,10,replace=False)
# # print(simplifier.nodes)
# # print(samples)
# x = x.simplify(samples=samples.tolist())
# x.dump_tables(nodes=nodes,edgesets=edgesets)
# # Based on Peter Ralph's example
# # at https://github.com/petrelharp/local_pca/blob/master/sims/msp/msp-add-mutation.py
# msp_rng = msprime.RandomGenerator(seed)
# mutations = msprime.MutationTable()
# sites = msprime.SiteTable()
# mutgen = msprime.MutationGenerator(msp_rng, theta/float(4*N)) # rho = theta
# mutgen.generate(nodes, edgesets, sites, mutations)
# x = msprime.load_tables(nodes=nodes,edgesets=edgesets,sites=sites,mutations=mutations)
# print(simplifier.nodes)
# print(nodes)
# for tree in x.trees():
#     for mut in tree.mutations():
#         print(tree.get_num_leaves(mut.node),tree.get_sample_size())
# print(sites.num_rows)

# np.random.seed(seed)
# seeds = np.random.randint(0,40000000,nreps)
# S = np.empty([len(seeds)],dtype=np.uint32)
# for i,ith_seed in enumerate(seeds):
#     ith_seed = 24331787
#     print("seed will be", ith_seed)
#     simplifier,atracker = bf.evolve_track_wrapper(popsize=N, rho=rho, seed = int(ith_seed), gc_interval=100)
# 
#     # Get a sample of size n = 10 
#     # This next set of steps is a little dumb:
#     # We are over-writing the data for our population,
#     # which is not desirable.
#     x = msprime.load_tables(nodes=simplifier.nodes, edgesets=simplifier.edgesets)
#     samples = np.random.randint(0,2*N,10,dtype=np.uint32)
#     x = x.simplify(samples=samples.tolist())
#     x.dump_tables(nodes=simplifier.nodes,edgesets=simplifier.edgesets)
#     # Based on Peter Ralph's example
#     # at https://github.com/petrelharp/local_pca/blob/master/sims/msp/msp-add-mutation.py
#     msp_rng = msprime.RandomGenerator(seed)
#     mutations = msprime.MutationTable()
#     sites = msprime.SiteTable()
#     mutgen = msprime.MutationGenerator(msp_rng, rho/float(4*N)) # rho = theta
#     mutgen.generate(simplifier.nodes, simplifier.edgesets, sites, mutations)
#     S[i] = sites.num_rows
#     print(S[i])
#     #print(sites)
# 
# print(S.mean())
