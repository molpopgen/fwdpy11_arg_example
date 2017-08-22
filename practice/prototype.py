# Prototype of our C++ implementation
# in Python using NumPy arrays and
# a VERY simple W-F simulations.
# Similar to what is in ftprime
# test suite, but with diploids.
# The data structures map to what
# I'm doing on the C++ side.
# NumPy arrays containing node_dy
# and edge_dt
# are built up the same way that
# I'm populating vector<node> and
# vector<edge>.

import numpy as np
import msprime
import time


node_dt = np.dtype([('id', np.uint32),
                    ('generation', np.float),
                    ('population', np.uint32)])

edge_dt = np.dtype([('left', np.float),
                    ('right', np.float),
                    ('parent', np.uint32),
                    ('child', np.uint32)])


def xover():
    breakpoint = np.random.sample()
    while breakpoint == 0.0 or breakpoint == 1.0:
        breakpoint = np.random.sample()
    return breakpoint


def wf(diploids, ngens):
    """
    For N diploids, the diploids list contains 2N values.
    For the i-th diploids, diploids[2*i] and diploids[2*i+1]
    represent the two chromosomes.

    We simulate constant N here, according to a Wright-Fisher
    scheme:

    1. Pick two parental indexes, each [0,N-1]
    2. The two chromosome ids are diploids[2*p] and diploids[2*p+1]
       for each parental index 'p'.
    3. We pass one parental index on to each offspring, according to
       Mendel, which means we swap parental chromosome ids 50% of the time.
    4. We do a single crossover for every mating, just so that there are
       a lot of breakpoints to handle
    """
    N = int(len(diploids) / 2)
    next_id = len(diploids)  # This will be the next unique ID to use
    assert(max(diploids) < next_id)
    # nodes = [Node(i, 0, 0) for i in diploids]  # Add nodes for ancestors
    nodes = np.array([(i, 0, 0) for i in diploids], dtype=node_dt)
    edges = np.empty([0], dtype=edge_dt)

    # We know there will be 2N new nodes added,
    # so we pre-allocate the space. We only need
    # to do this once b/c N is constant.
    tnodes = np.empty([2 * N], dtype=node_dt)

    tedges = np.empty([4 * N], dtype=edge_dt)
    for gen in range(ngens):
        # Empty offspring list.  We initialize
        # as a copy just to get the size right
        new_diploids = np.empty([len(diploids)], dtype=diploids.dtype)

        # Pick 2N parents:
        parents = np.random.randint(0, N, 2 * N)
        dip = int(0)
        edge_index = int(0)
        for parent1, parent2 in zip(parents[::2], parents[1::2]):
            # Pick two parents
            parents = np.random.randint(0, N, 2)
            # p1g1 = parent 1, gamete (chrom) 1, etc.:
            p1g1, p1g2 = diploids[2 * parent1], diploids[2 * parent1 + 1]
            p2g1, p2g2 = diploids[2 * parent2], diploids[2 * parent2 + 1]
            assert(p1g2 - p1g1 == 1)
            assert(p2g2 - p2g1 == 1)

            # Apply Mendel to p1g1 et al.
            mendel = np.random.random_sample(2)
            if mendel[0] < 0.5:
                p1g1, p1g2 = p1g2, p1g1
            if mendel[1] < 0.5:
                p2g1, p2g2 = p2g2, p2g1

            # We'll make every mating have 1 x-over
            # in each parent

            # Crossing-over and edges due to
            # contribution from parent 1
            breakpoint = xover()

            tedges[edge_index] = ((0.0, breakpoint, p1g1, next_id))
            tedges[edge_index + 1] = ((breakpoint, 1.0, p1g2, next_id))

            # Repeat process for parent 2's contribution
            breakpoint = xover()
            tedges[edge_index + 2] = ((0.0, breakpoint, p2g1, next_id + 1))
            tedges[edge_index + 3] = ((breakpoint, 1.0, p2g2, next_id + 1))

            # Add diploids
            new_diploids[2 * dip] = next_id
            new_diploids[2 * dip + 1] = next_id + 1

            # Add new nodes
            tnodes[2 * dip] = ((next_id, gen + 1, 0))
            tnodes[2 * dip + 1] = ((next_id + 1, gen + 1, 0))

            next_id += 2
            dip += 1
            edge_index += 4

        assert(edge_index == 4 * N)
        assert(len(new_diploids) == 2 * N)
        diploids = new_diploids
        assert(max(diploids) < next_id)
        nodes = np.concatenate([nodes, tnodes])
        edges = np.concatenate([edges, tedges])
    return (nodes, edges, diploids)


popsize = 100
diploids = np.array([i for i in range(2 * popsize)], dtype=np.uint32)
np.random.seed(42)
startsim = time.time()
ne = wf(diploids, 10 * popsize)
stopsim = time.time()

print("sim time took, ", stopsim - startsim)
nodes = ne[0]
edges = ne[1]


max_gen = max([i['generation'] for i in nodes])

samples = [i['id'] for i in nodes if i['generation'] == max_gen]

print(len(nodes), len(edges), len(samples))

print(nodes[:5])
print(nodes[-5:])
nodes['generation'] = nodes['generation'] - max_gen
nodes['generation'] = nodes['generation'] * -1.0
print(nodes[:5])
print(nodes[-5:])
# for i in nodes:
#     g = i['generation']
#     i.generation -= max_gen
#     i.generation *= -1

nt = msprime.NodeTable()
for i in nodes:
    flag = True if i['generation'] == 0 else False
    nt.add_row(flags=True, population=i['population'], time=i['generation'])

es = msprime.EdgesetTable()
for i in edges:
    es.add_row(left=i['left'], right=i['right'],
               parent=i['parent'], children=(i['child'],))

msprime.sort_tables(nodes=nt, edgesets=es)
x = msprime.load_tables(nodes=nt, edgesets=es)
x = x.simplify(samples=samples)

nt_s = msprime.NodeTable()
es_s = msprime.EdgesetTable()

x.dump_tables(nodes=nt_s, edgesets=es_s)

# print(nt)
# print(nt_s)
# print(es)
# print(es_s)
# print(samples)
