import fwdpy11_arg_example.evolve_without_simplify
import numpy as np
import msprime
x = fwdpy11_arg_example.evolve_without_simplify.evolve(1000, 10000)

sim_nodes = np.array(x.nodes, copy=False)
sim_edges = np.array(x.edges, copy=False)

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()

flags = np.ones(len(sim_nodes), dtype=np.uint32)
nodes.set_columns(
    flags=flags, population=sim_nodes['population'], time=sim_nodes['generation'])
edges.set_columns(left=sim_edges['left'], right=sim_edges['right'],
                  parent=sim_edges['parent'], child=sim_edges['child'])
msprime.sort_tables(nodes=nodes,edges=edges)

#Simplify w.r.to a sample of 100 haplotypes:
samples = np.random.choice(2000,100,replace=False)
msprime.simplify_tables(samples=samples.tolist(),nodes=nodes,edges=edges)
