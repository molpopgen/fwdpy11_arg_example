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
import sys

node_dt = np.dtype([('id', np.uint32),
                    ('generation', np.float),
                    ('population', np.int32)])

edge_dt = np.dtype([('left', np.float),
                    ('right', np.float),
                    ('parent', np.int32),
                    ('child', np.int32)])

mutation_dt = np.dtype([('position',np.float64), 
 						('node_id',np.int32)])

# Simulation with be popsize*SIMLEN generations
SIMLEN=20

class MockAncestryTracker(object):
    """
    Mimicking the public API of AncestryTracker.
    """
    __nodes = None
    __edges = None
    __mutations = None
    
    def __init__(self):
        self.nodes = np.empty([0], dtype=node_dt)
        self.edges = np.empty([0], dtype=edge_dt)
        self.mutations = np.empty([0], dtype=mutation_dt)
        
    @property
    def nodes(self):
        return self.__nodes

    @nodes.setter
    def nodes(self, value):
        self.__nodes = value

    @property
    def edges(self):
        return self.__edges

    @edges.setter
    def edges(self, value):
        self.__edges = value

    @property
    def mutations(self):
        return self.__mutations

    @mutations.setter
    def mutations(self, value):
        self.__mutations = value

def random_loc():
    breakpoint = np.random.sample()
    while breakpoint == 0.0 or breakpoint == 1.0:
        breakpoint = np.random.sample()
    return breakpoint


def wf(N, tracker, ngens):
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
       
    For each mating, each child chromosome gets 1 mutation in a random position
    """
    diploids = np.arange(2 * N, dtype=np.uint32)
    # so we pre-allocate the space. We only need
    # to do this once b/c N is constant.
    tracker.nodes = np.empty([2 * N * (ngens + 1)], dtype=node_dt)
    tracker.nodes['id'][:len(diploids)] = diploids
    tracker.nodes['generation'][:len(diploids)] = 0.0
    # Our simple WF sim makes 1 xover per parent
    # in each mating.  Thus, each offspring inherits
    # two new edges, and 4N new edges are formed each generation.
    tracker.edges = np.empty([ngens * 4 * N], dtype=edge_dt)
    edge_index = int(0)
    # Our simple WF sim makes 1 mutation per parent
    # in each mating.  Thus, each offspring inherits
    # two new mutations, and 4N new mutations are formed each generation.
    tracker.mutations = np.empty([ngens * 2 * N], dtype=mutation_dt)
    mutation_index = int(0)
    
    node_id = int(len(diploids))
    next_id = len(diploids)  # This will be the next unique ID to use
    assert(max(diploids) < next_id)
    for gen in range(ngens):
        # Empty offspring list.
        # We also use this to mark the "samples" for simplify
        new_diploids = np.empty([len(diploids)], dtype=diploids.dtype)

        # Pick 2N parents:
        parents = np.random.randint(0, N, 2 * N)
        assert(parents.max() < N)
        dip = int(0)  # dummy index for filling contents of new_diploids

        # Iterate over our chosen parents via fancy indexing.
        for parent1, parent2 in zip(parents[::2], parents[1::2]):
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
            breakpoint = random_loc()
            mutation_pos = random_loc()
            # Add the edges. Parent 1 will contribute [0,breakpoint)
            # from p1g1 and [breakpoint,1.0) from p1g2. Note that
            # we've already done the Mendel thing above. Both of these
            # segments are inherited by the next diploid, whose value
            # is next_id
            tracker.edges[edge_index] = (0.0, breakpoint, p1g1, next_id)
            tracker.edges[edge_index + 1] = (breakpoint, 1.0, p1g2, next_id)

            # Add mutations. Parent 1 will contribute 1 mutation at mutation_pos
            tracker.mutations[mutation_index] = (mutation_pos, next_id)
            # Update offspring container for
            # offspring dip, chrom 1:
            new_diploids[2 * dip] = next_id

            # Add new node for offpsring dip, chrom 1
            tracker.nodes[node_id] = (next_id, gen + 1, 0)

            # Repeat process for parent 2's contribution.
            # Stuff is now being inherited by node next_id + 1
            breakpoint = random_loc()
            mutation_pos = random_loc()
            
            tracker.edges[edge_index + 2] = (0.0, breakpoint, p2g1, next_id + 1)
            tracker.edges[edge_index + 3] = (breakpoint, 1.0, p2g2, next_id + 1)
            
            tracker.mutations[mutation_index+1] = (mutation_pos, next_id+1)

            new_diploids[2 * dip + 1] = next_id + 1

            tracker.nodes[node_id + 1] = (next_id + 1, gen + 1, 0)

            # python -O to turn this stuff off
            # if __debug__:
            #     print("generation ", gen, diploids.min(), diploids.max(),
            #           tracker.edges[edge_index:edge_index + 4])

            # Update our dummy variables.
            # We have two new nodes,
            # have used up two ids,
            # processed one diploid (offspring),
            # and added 4 edges
            node_id += 2
            next_id += 2
            mutation_index += 2
            dip += 1
            edge_index += 4

        assert(dip == N)
        assert(len(new_diploids) == 2 * N)
        diploids = new_diploids
        assert(max(diploids) < next_id)

    return (diploids)


def expensive_check(popsize, edges, nodes):
    """
    A brute-force post-hoc check of the
    nodes and edges that we generated
    in the simulation
    """
    assert(len(edges) == SIMLEN * popsize * 4 * popsize)

    # Check that all parent/child IDs are
    # in expected range.
    for gen in range(1, SIMLEN * popsize + 1):
        min_parental_id = 2 * popsize * (gen - 1)
        max_parental_id = 2 * popsize * (gen - 1) + 2 * popsize
        min_child_id = max_parental_id
        max_child_id = min_child_id + 2 * popsize

        node_gen_m1 = nodes[np.argwhere(
            nodes['generation'] == float(gen - 1)).flatten()]
        assert(len(node_gen_m1) == 2 * popsize)
        if any(i < min_parental_id or i >= max_parental_id for i in node_gen_m1['id']) is True:
            raise RuntimeError("generation", gen - 1, "confused")

        node_gen = nodes[np.argwhere(
            nodes['generation'] == float(gen)).flatten()]
        assert(len(node_gen) == 2 * popsize)
        if any(i < min_child_id or i >= max_child_id for i in node_gen['id']) is True:
            raise RuntimeError("generation", gen, "confused")

        edges_gen = edges[(gen - 1) * 4 * popsize:(gen - 1)
                          * 4 * popsize + 4 * popsize]

        if any(i not in node_gen_m1['id'] for i in edges_gen['parent']) is True:
            raise RuntimeError(
                "parent not found in expected slice of node table")

        if any(i not in node_gen['id'] for i in edges_gen['child']) is True:
            raise RuntimeError(
                "child not found in expected slice of node table")

        if any(i < min_parental_id or i >= max_parental_id for i in edges_gen['parent']) is True:
            raise RuntimeError("Bad parent")

        if any(i < min_child_id or i >= max_child_id for i in edges_gen['child']) is True:
            raise RuntimeError("Bad child")
    assert(float(gen) == nodes['generation'].max())


if __name__ == "__main__":
    popsize = int(sys.argv[1])
    theta = float(sys.argv[2])
    nsam = int(sys.argv[3])  # sample size to take and add mutations to
    seed = int(sys.argv[4])

    np.random.seed(seed)

    tracker = MockAncestryTracker()

    samples = wf(popsize, tracker, SIMLEN * popsize)

    # Check that our sample IDs are as expected:
    if __debug__:
        min_sample = SIMLEN * popsize * 2 * popsize
        max_sample = SIMLEN * popsize * 2 * popsize + 2 * popsize
        if any(i < min_sample or i >= max_sample for i in samples) is True:
            raise RuntimeError("Houston, we have a problem.")
    assert(np.array_equal(samples,np.arange(min_sample,max_sample,dtype=samples.dtype)))

    # Make local names for convenience
    nodes = tracker.nodes
    edges = tracker.edges
    mutas = tracker.mutations

    if __debug__:
        expensive_check(popsize, edges, nodes)

    max_gen = nodes['generation'].max()
    assert(int(max_gen) == SIMLEN * popsize)

    # Convert node times from forwards to backwards
    nodes['generation'] = nodes['generation'] - max_gen
    nodes['generation'] = nodes['generation'] * -1.0

    # Construct and populate msprime's tables
    flags = np.empty([len(nodes)], dtype=np.uint32)
    flags.fill(1)
    nt = msprime.NodeTable()
    nt.set_columns(flags=flags,
                   population=nodes['population'],
                   time=nodes['generation'])

    es = msprime.EdgeTable()
    es.set_columns(left=edges['left'],
                   right=edges['right'],
                   parent=edges['parent'],
                   child=edges['child'])
    
    st = msprime.SiteTable()
    st.set_columns(position=mutas['position'],
                   ancestral_state=np.zeros(len(mutas['position']),np.int8),
                   ancestral_state_length=np.ones(len(mutas['position']),np.uint32))
    
    mt = msprime.MutationTable()
    mt.set_columns(site=np.arange(len(mutas['node_id']),dtype=np.int32),
                   node=mutas['node_id'],
                   derived_state=np.ones(len(mutas['node_id']),np.int8),
                   derived_state_length=np.ones(len(mutas['node_id']),np.uint32))

    # Sort
    msprime.sort_tables(nodes=nt, edges=es, sites=st, mutations=mt)
    print("num total mutations: ", st.num_rows)
     
    # Simplify: this is where the magic happens
    ## PLR: since these tables aren't valid, you gotta use simplify_tables, not load them into a tree sequence
    nt_c = nt.copy()
    es_c = es.copy()
    st_c = st.copy()
    mt_c = mt.copy()
    msprime.simplify_tables(samples=samples.tolist(), nodes=nt_c, edges=es_c, sites=st_c, mutations=mt_c)
    print("num simplified mutations: ", st_c.num_rows)
    # Create a tree sequence
    x = msprime.load_tables(nodes=nt_c, edges=es_c, sites=st_c, mutations=mt_c)   
    
    print(max(mt_c.node))
    print(nt_c.num_rows)
    
    nt_s = nt_c.copy()
    es_s = es_c.copy()
    st_s = st_c.copy()
    mt_s = mt_c.copy()
       
    nsam_samples = np.random.choice(2 * popsize, nsam, replace=False)
    xs = x.simplify(nsam_samples.tolist())
    xs.dump_tables(nodes=nt_s, edges=es_s, sites=st_s, mutations=mt_s)
    
    print(max(mt_s.node))
    print(nt_s.num_rows)
