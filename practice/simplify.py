#mimic more closely what will come out of 
#the C++ back-end
import msprime
import numpy as np


n = msprime.NodeTable()
sv = [True, True, True, True, True, True, True]
tv = [0.0, 0.0, 0.0, 0.4, 0.5, 0.7, 1.0]
pv = [0, 0, 0, 0, 0, 0, 0]
n = msprime.NodeTable()
n.set_columns(flags=sv, population=pv, time=tv)
print(n)

left = [0.2, 0.2, 0.0, 0.0, 0.2, 0.2, 0.8, 0.8, 0.8, 0.8, 0.0, 0.0]
right = [0.8, 0.8, 0.2, 0.2, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 0.2, 0.2]
parent = [3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 6, 6]
# children = [(0,2),(1,2),(1,3),(0,4),(0,4)]
children = [0, 2, 1, 2, 1, 3, 1, 2, 0, 4, 0, 4]

e = msprime.EdgesetTable()

for l, r, p, c in zip(left, right, parent, children):
    e.add_row(left=l, right=r, parent=p, children=(c,))

print(e)
msprime.sort_tables(nodes=n,edgesets=e)
x = msprime.load_tables(nodes=n,edgesets=e)
x = x.simplify(samples=[0,1,2])
x.dump_tables(nodes=n,edgesets=e)
print(n)
print(e)
#make some fake nodes 
