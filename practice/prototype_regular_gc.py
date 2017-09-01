# Advanced prototype:
# Does GC at user-specified intervals
# Random number of edges due to modeling
# recombination as a Poisson process.

import numpy as np
import msprime
import sys
import argparse

node_dt = np.dtype([('id', np.uint32),
                    ('generation', np.float),
                    ('population', np.int32)])

edge_dt = np.dtype([('left', np.float),
                    ('right', np.float),
                    ('parent', np.int32),
                    ('child', np.int32)])

# Simulation with be popsize*SIMLEN generations
SIMLEN = 20


class MockAncestryTracker(object):
    """
    Mimicking the public API of AncestryTracker.
    """
    __nodes = None
    __edges = None
    __samples = None

    def __init__(self):
        self.nodes = np.empty([0], dtype=node_dt)
        self.edges = np.empty([0], dtype=edge_dt)
        self.samples = None

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
    def samples(self):
        return self.__samples

    @samples.setter
    def samples(self, value):
        self.__samples = value

    def update_data(self, new_nodes, new_edges, new_samples):
        self.nodes = np.insert(self.nodes, len(self.nodes), new_nodes)
        self.edges = np.insert(self.edges, len(self.edges), new_edges)
        self.samples = new_samples

    def convert_time(self):
        self.nodes['generation'] -= self.nodes['generation'].max()
        self.nodes['generation'] *= -1.0

    def reset_data(self):
        """
        Call this after garbage collection.
        Any necessary post-processing is done here.
        """
        self.nodes = np.empty([0], dtype=node_dt)
        self.edges = np.empty([0], dtype=edge_dt)

    def post_gc_cleanup(self, gc_rv):
        if gc_rv[0] is True:
            self.reset_data()


class ARGsimplifier(object):
    """
    Mimicking the API if our Python
    class to collect simulated
    results and process them via msprime
    """
    __nodes = None
    __edges = None
    __gc_interval = None
    __last_gc_time = None

    def __init__(self, gc_interval=None):
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgesetTable()
        self.gc_interval = gc_interval
        self.last_gc_time = 0

    def simplify(self, generation, tracker):
        """
        Details of taking new data, appending, and
        simplifying.

        :return: length of simplifed node table, which is next_id to use
        """
        # Update time in current nodes.
        # Is this most effficient method?
        dt = generation - self.last_gc_time
        self.nodes.set_columns(flags=self.nodes.flags,
                               population=self.nodes.population,
                               time=self.nodes.time + dt)

        # Create "flags" for new nodes.
        # This is much faster than making a list
        flags = np.empty([len(tracker.nodes)], dtype=np.uint32)
        flags.fill(1)

        tracker.convert_time()

        self.nodes.append_columns(flags=flags,
                                  population=tracker.nodes['population'],
                                  time=tracker.nodes['generation'])
        self.edges.append_columns(left=tracker.edges['left'],
                                  right=tracker.edges['right'],
                                  parent=tracker.edges['parent'],
                                  children=tracker.edges['child'],
                                  children_length=[1] * len(tracker.edges))

        msprime.sort_tables(nodes=self.nodes, edgesets=self.edges)
        msprime.simplify_tables(samples=tracker.samples.tolist(),
                                nodes=self.nodes, edgesets=self.edges)
        return self.nodes.num_rows

    def __call__(self, generation, tracker):
        if generation > 0 and generation % self.gc_interval == 0:
            next_id = self.simplify(generation, tracker)
            self.last_gc_time = generation
            return (True, next_id)

        return (False, None)

    @property
    def gc_interval(self):
        return self.__gc_interval

    @gc_interval.setter
    def gc_interval(self, value):
        self.__gc_interval = int(value)

    @property
    def last_gc_time(self):
        return self.__last_gc_time

    @last_gc_time.setter
    def last_gc_time(self, value):
        self.__last_gc_time = int(value)

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


def xover(rate):
    nbreaks = np.random.poisson(rate)
    if nbreaks == 0:
        return np.empty([0], dtype=np.float)
    rv = np.random.random_sample(nbreaks)
    # np.unique both sorts and purges
    # duplicate values. (Duplicate values
    # are functionally a double x-over)
    rv = np.unique(rv)
    # We "cap" the set of breakpoints
    # as in fwdpp:
    rv = np.insert(rv, len(rv), np.finfo(np.float).max)
    return rv


def split_breakpoints(breakpoints):
    """
    Take the breakpoints from a meiosis,
    and return them as segments contributed
    by gamete 1 and gamete 2

    Note: bug source could be here, if breakpoints[0] == 0.0
    """
    s1 = np.array([(0.0, breakpoints[0])], dtype=[
                  ('left', np.float), ('right', np.float)])
    s2 = np.empty([0], dtype=s1.dtype)
    for i in range(1, len(breakpoints)):
        a = breakpoints[i - 1]
        b = breakpoints[i] if i < len(breakpoints) - 1 else 1.0
        assert(a != b)
        if i % 2 == 0.:
            s1 = np.insert(s1, len(s1), (a, b))
        else:
            s2 = np.insert(s2, len(s2), (a, b))
    return (s1, s2)


def handle_recombination_update(offspring_index, parental_id1,
                                parental_id2, edges, breakpoints):
    if len(breakpoints) == 0:
        edges.append((0.0, 1.0, parental_id1, offspring_index))
        # edges = np.insert(edges, len(edges),
        #                  (0.0, 1.0, parental_id1, offspring_index))
        return edges

    split = split_breakpoints(breakpoints)
    for i, j in split[0]:
        edges.append((i, j, parental_id1, offspring_index))
        # edges = np.insert(edges, len(edges),
        #                   (i, j, parental_id1, offspring_index))
    for i, j in split[1]:
        edges.append((i, j, parental_id2, offspring_index))
        # edges = np.insert(edges, len(edges),
        #                   (i, j, parental_id2, offspring_index))
    return edges


def wf(N, simplifier, tracker, recrate, ngens):
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
    4. Crossing over is a Poisson process, and the code used here (above)
       is modeled after that in our C++ implementation.
    """
    diploids = np.arange(2 * N, dtype=np.uint32)
    # so we pre-allocate the space. We only need
    # to do this once b/c N is constant.
    tracker.nodes = np.array([(i, 0, 0) for i in diploids], dtype=node_dt)

    next_id = len(diploids)  # This will be the next unique ID to use
    assert(max(diploids) < next_id)
    for gen in range(ngens):
        # Let's see if we will do some GC:
        gc_rv = simplifier(gen, tracker)
        # If so, let the tracker clean up:
        tracker.post_gc_cleanup(gc_rv)
        if gc_rv[0] is True:
            assert(len(tracker.nodes)==0)
            assert(len(tracker.edges)==0)
            # If we did GC, we need to reset
            # some variables.  Internally,
            # when msprime simplifies tables,
            # the 2N tips are entries 0 to 2*N-1
            # in the NodeTable, hence the 
            # re-assignment if diploids.
            # We can also reset next_id to the
            # length of the current NodeTable,
            # keeping risk of integer overflow 
            # to a minimum.
            next_id = int(gc_rv[1])
            diploids = np.arange(2 * N, dtype=np.uint32)

        # Empty offspring list.
        # We also use this to mark the "samples" for simplify
        new_diploids = np.empty([len(diploids)], dtype=diploids.dtype)

        # Store temp nodes for this generation.
        # b/c we add 2N nodes each geneartion,
        # We can allocate it all at once.
        nodes = np.empty([2 * N], dtype=node_dt)
        nodes['population'].fill(0)
        nodes['generation'] = gen + 1
        nodes['id'] = np.arange(
            start=next_id, stop=next_id + 2 * N, dtype=nodes['id'].dtype)

        # Store temp edges for this generation
        edges = []  # np.empty([0], dtype=edge_dt)
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

            breakpoints = xover(recrate)

            edges = handle_recombination_update(
                next_id, p1g1, p1g2, edges, breakpoints)

            # Update offspring container for
            # offspring dip, chrom 1:
            new_diploids[2 * dip] = next_id

            # Repeat process for parent 2's contribution.
            # Stuff is now being inherited by node next_id + 1
            breakpoints = xover(recrate)
            edges = handle_recombination_update(
                next_id + 1, p2g1, p2g2, edges, breakpoints)

            new_diploids[2 * dip + 1] = next_id + 1

            # Update our dummy variables.
            next_id += 2
            dip += 1

        assert(dip == N)
        assert(len(new_diploids) == 2 * N)
        diploids = new_diploids
        tracker.update_data(nodes, edges, diploids)
        assert(max(diploids) < next_id)

    return (diploids)


def parse_args():
    dstring = "Prototype implementation of ARG tracking and regular garbage collection."
    parser = argparse.ArgumentParser(description=dstring,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--popsize', '-N', type=int,
                        default=500, help="Diploid population size")
    parser.add_argument('--theta', '-T', type=float, default=10.0, help="4Nu")
    parser.add_argument('--rho', '-R', type=float, default=10.0, help="4Nr")
    parser.add_argument('--nsam', '-n', type=int, default=10,
                        help="Sample size (in chromosomes).")
    parser.add_argument('--seed', '-S', type=int, default=42, help="RNG seed")
    parser.add_argument('--gc', '-G', type=int,
                        default=100, help="GC interval")

    return parser


if __name__ == "__main__":
    parser = parse_args()
    args = parser.parse_args(sys.argv[1:])

    np.random.seed(args.seed)

    simplifier = ARGsimplifier(args.gc)
    tracker = MockAncestryTracker()
    recrate = args.rho / float(4 * args.popsize)
    samples = wf(args.popsize, simplifier, tracker,
                 recrate, SIMLEN * args.popsize)

    if len(tracker.nodes) > 0:  # Then there's stuff that didn't get GC'd
        simplifier.simplify(SIMLEN*args.popsize, tracker)

    # Local names for convenience.
    # I copy the tables here, too,
    # because I think that will be
    # done in practice: you will
    # often want to simplify and
    # ARG down to a smaller sample
    # but still have the complete
    # history of the pop'n.
    nodes = simplifier.nodes.copy()
    edges = simplifier.edges.copy()

    nsam_samples = np.random.choice(2 * args.popsize, args.nsam, replace=False)
    msprime.simplify_tables(samples=nsam_samples.tolist(),
                            nodes=nodes, edgesets=edges)
    msp_rng = msprime.RandomGenerator(args.seed)
    mutations = msprime.MutationTable()
    sites = msprime.SiteTable()
    mutgen = msprime.MutationGenerator(
        msp_rng, args.theta / float(4 * args.popsize))
    mutgen.generate(nodes, edges, sites, mutations)
    print(sites.num_rows)
