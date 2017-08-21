# Messing around with msprime's
# Python API for nodes and edges
import msprime

n = msprime.NodeTable()
sv = [True, True, True, False, False, False, False]
tv = [0.0, 0.0, 0.0, 0.4, 0.5, 0.7, 1.0]
pv = [0, 0, 0, 0, 0, 0, 0]
for s, t, p in zip(sv, tv, pv):
    n.add_row(flags=s, population=p, time=t)


print(n)

n = msprime.NodeTable()
n.set_columns(flags=sv, population=pv, time=tv)

print(n)

left = [0.2, 0.0, 0.2, 0.8, 0.8, 0.0]
right = [0.8, 0.2, 0.8, 1.0, 1.0, 0.2]
parent = [3, 4, 4, 4, 5, 6]
# children = [(0,2),(1,2),(1,3),(0,4),(0,4)]
children = [(0, 2), (1, 2), (1, 3), (1, 2), (0, 4), (0, 4)]

e = msprime.EdgesetTable()

# e.set_columns(left=left,right=right,parent=parent,children=children,children_length=len(children))
for l, r, p, c in zip(left, right, parent, children):
    e.add_row(left=l, right=r, parent=p, children=c)


print(e)

msprime.sort_tables(nodes=n, edgesets=e)

print(n)

ts = msprime.load_tables(nodes=n, edgesets=e)

print(ts)

print(ts.samples())
