# Messing around with msprime's
# Python API for nodes and edges
import msprime

n = msprime.NodeTable()
sv = [True, True, True, False, False, False, False]
tv = [0.0, 0.0, 0.0, 0.4, 0.5, 0.7, 1.0]
pv = [0, 0, 0, 0, 0, 0, 0]
for s, t, p in zip(sv, tv, pv):
    n.add_row(flags=s, population=p, time=t)

n = msprime.NodeTable()
n.set_columns(flags=sv, population=pv, time=tv)

left = [0.2, 0.0, 0.2, 0.8, 0.8, 0.0]
right = [0.8, 0.2, 0.8, 1.0, 1.0, 0.2]
parent = [3, 4, 4, 4, 5, 6]
children = [0, 2, 1, 2, 1, 3, 1, 2, 0, 4, 0, 4]

#We swap some data around, but the tree is actuall the same.
left2 = [0.2, 0.0, 0.8, 0.2, 0.8, 0.0]
right2 = [0.8, 0.2, 1.0, 0.8, 1.0, 0.2]
children2 = [0, 2,
             1, 2,
             1, 2,
             1, 3,
             0, 4,
             0, 4]


def make_EdgesetTable(left, right, parent, children):
    e = msprime.EdgesetTable()
    e.set_columns(left=left, right=right, parent=parent,
                  children=children, children_length=[2] * len(parent))
    return e


def make_tree_add_mutations(nodes, edges, mutrate):
    rng = msprime.RandomGenerator(42)
    m = msprime.MutationTable()
    s = msprime.SiteTable()
    mg = msprime.MutationGenerator(rng, mutrate)
    mg.generate(nodes, edges, s, m)
    rv = msprime.load_tables(nodes=nodes, edgesets=edges, sites=s, mutations=m)
    return (rv, s)


e = make_EdgesetTable(left, right, parent, children)
e2 = make_EdgesetTable(left2, right2, parent, children2)

print(e)
print(e2)

m = make_tree_add_mutations(n, e, 10)
m2 = make_tree_add_mutations(n, e2, 10)

print("first:")
print(m[1])
print("second:")
print(m2[1])
