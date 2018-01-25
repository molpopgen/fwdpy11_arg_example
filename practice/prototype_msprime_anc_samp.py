import msprime
import numpy as np

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()
mutations = msprime.MutationTable()
sites = msprime.SiteTable()

prior_ts = msprime.simulate(sample_size=2 * 500, Ne=2 *
                            500, mutation_rate=10 / float(4 * 500), random_seed=10)
prior_ts.dump_tables(nodes=nodes, edges=edges,
                     sites=sites, mutations=mutations)

flags = np.empty([len(nodes)], dtype=np.uint32)
flags.fill(1)

nodes.set_columns(flags=flags, population=nodes.population, time=nodes.time)
samples = [1992, 1989, 1995] + np.arange(100).tolist()
node_map = msprime.simplify_tables(
    samples=samples, nodes=nodes, edges=edges, sites=sites, mutations=mutations)

print(mutations)

x = msprime.load_tables(nodes=nodes, edges=edges,
                        sites=sites, mutations=mutations)
t = next(x.trees())
t.draw("../../test_tree_1.svg", width=1500, height=1500)

samples = [0, 1, 2] + (np.arange(10) + 3).tolist()

flags = np.empty([len(nodes)], dtype=np.uint32)
flags.fill(1)
nodes.set_columns(flags=flags, population=nodes.population, time=nodes.time)

node_map = msprime.simplify_tables(
    samples=samples, nodes=nodes, edges=edges, sites=sites, mutations=mutations)
print(mutations)


x = msprime.load_tables(nodes=nodes, edges=edges,
                        sites=sites, mutations=mutations)
t = next(x.trees())
t.draw("../../test_tree_2.svg", width=1500, height=1500)
