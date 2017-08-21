# A VERY simple W-F simulations.
# Similar to what is in ftprime
# test suite, but with diploids.
# The data structures map to what
# I'm doing on the C++ side.
# A list[Node] and list[Edge]
# are built up the same way that
# I'm populating vector<node> and
# vector<edge>.

import numpy as np
import msprime


class Node(object):
    """
    Tree nodes
    """
    id = None
    generation = None
    population = None

    def __init__(self, id, generation, population):
        self.id = id
        self.generation = generation
        self.population = population

    def __repr__(self):
        rv = "Node("
        rv += str(self.id) + ','
        rv += str(self.generation) + ','
        rv += str(self.population) + ')'
        return rv


class Edge(object):
    """
    Simple representation of edges.
    We expect msprime to "collect"
    children during the simplification
    steps.
    """
    left = None
    right = None
    parent = None
    child = None

    def __init__(self, left, right, parent, child=None):
        self.left = left
        self.right = right
        self.parent = parent
        self.child = child

    def __repr__(self):
        rv = "Edge("
        rv += str(self.left) + ','
        rv += str(self.right) + ','
        rv += str(self.parent) + ','
        rv += str(self.child) + ')'
        return rv


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
    nodes = [Node(i, 0, 0) for i in diploids]  # Add nodes for ancestors
    edges = []
    for gen in range(ngens):
        # Empty offspring list.  We initialize
        # as a copy just to get the size right
        new_diploids = np.array(diploids, copy=True)
        for dip in range(N):
            # Pick two parents
            parents = np.random.randint(0, N, 2)
            # p1g1 = parent 1, gamete (chrom) 1, etc.:
            p1g1, p1g2 = diploids[2 * parents[0]], diploids[2 * parents[0] + 1]
            p2g1, p2g2 = diploids[2 * parents[1]], diploids[2 * parents[1] + 1]
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

            # Crossing-over and Edges due to
            # contribution from parent 1
            breakpoint = xover()

            edges.append(Edge(0.0, breakpoint, p1g1, next_id))
            edges.append(Edge(breakpoint, 1.0, p1g2, next_id))

            # Repeat process for parent 2's contribution
            breakpoint = xover()
            edges.append(Edge(0.0, breakpoint, p2g1, next_id + 1))
            edges.append(Edge(breakpoint, 1.0, p2g2, next_id + 1))

            new_diploids[2*dip] = next_id
            new_diploids[2*dip+1] = next_id+1
            next_id += 2

        assert(len(new_diploids) == 2 * N)
        diploids = new_diploids
        assert(max(diploids) < next_id)
        for i in diploids:
            nodes.append(Node(i, gen + 1, 0))

    return (nodes, edges, diploids)


popsize = 100
diploids = np.array([i for i in range(2 * popsize)], dtype=np.uint32)
np.random.seed(42)
ne = wf(diploids, 10*popsize)

nodes = ne[0]
edges = ne[1]

max_gen = max([i.generation for i in nodes])

samples = [i.id for i in nodes if i.generation == max_gen]

print(len(nodes), len(edges), len(samples))

for i in nodes:
    g = i.generation
    i.generation -= max_gen
    i.generation *= -1

nt = msprime.NodeTable()
for i in nodes:
    flag = True if i.generation == 0 else False
    nt.add_row(flags=True, population=i.population, time=i.generation)

es = msprime.EdgesetTable()
for i in edges:
    es.add_row(left=i.left, right=i.right,
               parent=i.parent, children=(i.child,))

# The tables are already sorted...
msprime.sort_tables(nodes=nt, edgesets=es)
x = msprime.load_tables(nodes=nt, edgesets=es)
x = x.simplify(samples=samples)

nt_s = msprime.NodeTable()
es_s = msprime.EdgesetTable()

x.dump_tables(nodes=nt_s, edgesets=es_s)

# print(nt)
print(nt_s)
# print(es)
# print(es_s)
# print(samples)
