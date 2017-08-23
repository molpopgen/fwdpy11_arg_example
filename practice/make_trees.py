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

e = msprime.EdgesetTable()
print("Jerome's example:")
left = [0, 1]
right = [1, 2]
parent = [10, 11]
children_length = [2, 3]
children = [1, 2, 3, 4, 5]
e.set_columns(left=left, right=right, parent=parent, children=children_length, children_length=children_length)
print(e)


left = [0.2, 0.0, 0.2, 0.8, 0.8, 0.0]
right = [0.8, 0.2, 0.8, 1.0, 1.0, 0.2]
parent = [3, 4, 4, 4, 5, 6]
# children = [(0,2),(1,2),(1,3),(0,4),(0,4)]
children = [0, 2, 1, 2, 1, 3, 1, 2, 0, 4, 0, 4]
#children = [(0, 2), (1, 2), (1, 3), (1, 2), (0, 4), (0, 4)]

print([2]*len(children))
e.set_columns(left=left,right=right,parent=parent,children=children,children_length=[2]*len(parent))
# for l, r, p, c in zip(left, right, parent, children):
#     e.add_row(left=l, right=r, parent=p, children=c)


print(e)



msprime.sort_tables(nodes=n, edgesets=e)

print(n)

ts = msprime.load_tables(nodes=n, edgesets=e)

print(ts)

print(ts.samples())

# Let's make an example that looks more like what comes out of 
# forward sims

left = [0.0, 0.5, 0.0, 0.2]
right = [0.5, 1.0 , 0.2, 1.0]
parent = [0, 0, 1, 1]
children = [2, 2, 3, 3]
children_length = [1]*4

e.set_columns(left=left,right=right,parent=parent,children=children,children_length=children_length)
print(e)
