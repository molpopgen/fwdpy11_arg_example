# A VERY simple W-F simulations.
# Similar to what is in ftprime
# test suite, but with diploids

import numpy as np
import msprime


class Node(object):
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


def wf(diploids, ngens):
    N = int(len(diploids) / 2)
    next_id = len(diploids)
    assert(max(diploids) < next_id)
    nodes = [Node(i, 0, 0) for i in diploids]
    edges = []
    for gen in range(ngens):
        new_diploids = []
        for dip in range(N):
            parents = np.random.randint(0, N, 2)
            p1g1, p1g2 = diploids[2 * parents[0]], diploids[2 * parents[0] + 1]
            assert(p1g2-p1g1 == 1)
            p2g1, p2g2 = diploids[2 * parents[1]], diploids[2 * parents[1] + 1]
            assert(p2g2-p2g1 == 1)
            mendel = np.random.random_sample(2)
            switch1, switch2 = 0, 0
            if mendel[0] < 0.5:
                switch1 = 1
            if mendel[1] < 0.5:
                switch2 = 1

            p1_chrom_label = p1g1 if switch1 == 0 else p1g2
            p2_chrom_label = p2g1 if switch2 == 0 else p2g2

            # We'll make every mating have 1 x-over
            breakpoint = np.random.random_sample()
            while breakpoint == 0.0 or breakpoint == 1.0:
                breakpoint = np.random.random_sample()

            edges.append(Edge(0.0, breakpoint, p1_chrom_label, next_id))
            edges.append(Edge(breakpoint, 1.0, p2_chrom_label, next_id + 1))
            new_diploids.append(next_id)
            new_diploids.append(next_id + 1)
            next_id += 2

        assert(len(new_diploids) == 2*N)
        diploids = new_diploids
        assert(max(diploids) < next_id)
        for i in diploids:
            nodes.append(Node(i, gen + 1, 0))

    return (nodes, edges, diploids)


diploids = []
for i in range(10):
    diploids.append(2 * i)
    diploids.append(2 * i + 1)
np.random.seed(42)
ne = wf(diploids, 10)

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
print(es_s)
# print(samples)
