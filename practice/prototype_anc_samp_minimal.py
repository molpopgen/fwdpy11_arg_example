# Minimal prototype:
# Does GC at user-specified intervals, ancestral sampling, prior history

import numpy as np
import msprime
import sys
import argparse

# Simulation with be popsize*SIMLEN generations
SIMLEN = 20

class MockAncestryTracker(object):
    samples = None
    anc_samples = None

    def __init__(self):
        self.samples = np.empty([0], dtype=np.uint32)
        self.anc_samples = np.empty([0], dtype=np.uint32)

    def update_samples(self, new_samples, new_anc_samples):
        """
        new_samples is the list of node IDs corresponding to the current children.
        """
        self.samples = new_samples
        self.anc_samples = np.insert(self.anc_samples, len(self.anc_samples), new_anc_samples)

class ARGsimplifier(object):
    """
    Mimicking the API of our Python
    class to collect simulated
    results and process them via msprime
    """
    nodes = None
    edges = None
    gc_interval = None

    def __init__(self, gc_interval=None):
        self.nodes = msprime.NodeTable()
        self.edges = msprime.EdgeTable()
        self.gc_interval = gc_interval

    def simplify(self, generation, ngens, tracker):
        """
        Details of simplifying.

        :return: length of simplifed node table, which is next_id to use
        """
        if(generation == 0):
            prior_ts = msprime.simulate(sample_size=2 * args.popsize, Ne=2 * args.popsize, random_seed=args.seed)
            prior_ts.dump_tables(nodes=self.nodes, edges=self.edges)
            self.nodes.set_columns(flags=self.nodes.flags, population=self.nodes.population, time=self.nodes.time + ngens)
            
            return self.nodes.num_rows

        # Sort and simplify
        msprime.sort_tables(nodes=self.nodes, edges=self.edges)
        # set guards against duplicate node ids when an ancestral sample generation overlaps with a gc generation
        # sorting the set in reverse order ensures that generational samples occur *before* ancestral samples in the node table,
        # making bookkeeping easier during the WF (parents of the next generation are guaranteed to be 0-2N in the node table)
        all_samples = sorted(set(tracker.anc_samples.tolist() + tracker.samples.tolist()), reverse=True)
        node_map = msprime.simplify_tables(samples=all_samples, nodes=self.nodes, edges=self.edges)

        tracker.anc_samples = np.array([node_map[int(node_id)] for node_id in tracker.anc_samples])
        # Return length of NodeTable,
        # which can be used as next offspring ID
        return self.nodes.num_rows

def wf(N, simplifier, tracker, anc_sample_gen, ngens):
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
    """
    next_id = simplifier.simplify(0, ngens, tracker)  #  Based on prior history, this will be the next unique ID to use
    diploids = np.arange(2 * N, dtype=np.uint32)
    
    ancestral_gen_counter = 0
    for gen in range(ngens):
        ancestral_samples = []

        # Empty offspring list.
        # We also use this to mark the "samples" for simplify
        new_diploids = np.empty([len(diploids)], dtype=diploids.dtype)

        if(ancestral_gen_counter < len(anc_sample_gen) and gen + 1 == anc_sample_gen[ancestral_gen_counter][0]):
            ran_samples = np.random.choice(int(N), int(anc_sample_gen[ancestral_gen_counter][1]), replace=False)
            # while sorting to get diploid chromosomes next to each other isn't strictly necessary,
            # they will be sorted (in reverse order) before simplication anyway, no need to do it here
            ancestral_samples = np.concatenate((2*ran_samples + next_id, 2*ran_samples + 1 + next_id))
            ancestral_gen_counter += 1

        # Store nodes for this generation.
        for i in range(2*N): simplifier.nodes.add_row(1,(ngens - (gen+1)),0)

        # Pick 2N parents:
        parents = np.random.randint(0, N, 2 * N)
        dip = int(0)  # dummy index for filling contents of new_diploids

        # Iterate over our chosen parents via fancy indexing.
        for parent1, parent2 in zip(parents[::2], parents[1::2]):
            mendel = np.random.random_sample(2)
            p1 = diploids[2 * parent1 + (mendel[0] < 0.5)]
            p2 = diploids[2 * parent2 + (mendel[1] < 0.5)]

            simplifier.edges.add_row(0.0, 1.0, p1, next_id)
            # Update offspring container for
            # offspring dip, chrom 1:
            new_diploids[2 * dip] = next_id

            # Repeat process for parent 2's contribution.
            # Stuff is now being inherited by node next_id + 1
            simplifier.edges.add_row(0.0, 1.0, p2, next_id+1)

            new_diploids[2 * dip + 1] = next_id + 1

            # Update our dummy variables.
            next_id += 2
            dip += 1

        diploids = new_diploids
        tracker.update_samples(diploids, ancestral_samples)
        # Let's see if we will do some GC:
        if ((gen+1) % simplifier.gc_interval == 0) or (gen == ngens):
            # If we do GC, we need to reset
            # some variables.  Internally,
            # when msprime simplifies tables,
            # the 2N tips are entries 0 to 2*N-1
            # in the NodeTable, hence the
            # re-assignment if diploids.
            # We can also reset next_id to the
            # length of the current NodeTable,
            # keeping risk of integer overflow
            # to a minimum.
            next_id = simplifier.simplify(gen, ngens, tracker)
            diploids = np.arange(2 * N, dtype=np.uint32)

    return (diploids)
    
def parse_args():
    dstring = "Prototype implementation of ARG tracking and regular garbage collection."
    parser = argparse.ArgumentParser(description=dstring, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--popsize', '-N', type=int, default=500, help="Diploid population size")
    parser.add_argument('--nsam', '-n', type=int, default=5, help="Sample size (in diploids).")
    parser.add_argument('--seed', '-S', type=int, default=42, help="RNG seed")
    parser.add_argument('--gc', '-G', type=int, default=100, help="GC interval")
    
    return parser

if __name__ == "__main__":
    parser = parse_args()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)

    simplifier = ARGsimplifier(args.gc)
    tracker = MockAncestryTracker()
    ngens = SIMLEN * args.popsize
    anc_sample_gen = [(ngens * (i + 1) / SIMLEN, max(round(args.popsize / 200), 1)) for i in range(SIMLEN - 2)]
    samples = wf(args.popsize, simplifier, tracker, anc_sample_gen, ngens)

    ran_samples = np.random.choice(args.popsize, args.nsam, replace=False)
    nsam_samples = sorted(np.concatenate((2*ran_samples, 2*ran_samples + 1)))
    all_samples = nsam_samples + tracker.anc_samples.tolist()
    node_map = msprime.simplify_tables(samples=all_samples, nodes=simplifier.nodes, edges=simplifier.edges)