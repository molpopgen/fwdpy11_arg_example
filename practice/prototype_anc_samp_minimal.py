# Minimal prototype:
# Does GC at user-specified intervals, ancestral sampling, prior history

import numpy as np
import msprime
import sys
import argparse

# Simulation with be popsize*SIMLEN generations
SIMLEN = 20

class ARGsimplifier(object):
    """
    Mimicking the API of our Python
    class to collect simulated
    results and process them via msprime
    """
    def __init__(self, ngens, gc_interval=None):
        self.nodes, self.edges = msprime.NodeTable(), msprime.EdgeTable()
        prior_ts = msprime.simulate(sample_size=2 * args.popsize, Ne=2 * args.popsize, random_seed=args.seed)
        prior_ts.dump_tables(nodes=self.nodes, edges=self.edges)
        self.nodes.set_columns(flags=self.nodes.flags, population=self.nodes.population, time=self.nodes.time + ngens)
        self.gc_interval = gc_interval

    def simplify(self, generation, ngens, samples, anc_samples):
        """
        Details of simplifying.

        :return: length of simplifed node table, which is next_id to use
        """
        diploids = np.arange(2 * args.popsize, dtype=np.uint32)
        if(generation > 0): 
            # Sort and simplify
            msprime.sort_tables(nodes=self.nodes, edges=self.edges)
            # set guards against duplicate node ids when an ancestral sample generation overlaps with a gc generation
            # sorting the set in reverse order ensures that generational samples occur *before* ancestral samples in the node table,
            # making bookkeeping easier during the WF (parents of the next generation are guaranteed to be 0-2N in the node table)
            all_samples = sorted(set(anc_samples.tolist() + samples.tolist()), reverse=True)
            node_map = msprime.simplify_tables(samples=all_samples, nodes=self.nodes, edges=self.edges)
            anc_samples = np.array([node_map[int(node_id)] for node_id in anc_samples])
            # Return length of NodeTable, which can be used as next offspring ID and reset diploids (parent gen to 0, 2N)
        return self.nodes.num_rows, diploids, diploids + self.nodes.num_rows, anc_samples

def wf(N, simplifier, anc_sample_gen, ngens):
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
    next_id, diploids, new_diploids, ancestral_samples = simplifier.simplify(0, ngens, np.empty([0], dtype=np.uint32), np.empty([0], dtype=np.uint32))  #  Based on prior history, this will be the next unique ID to use
    ancestral_gen_counter = 0
    for gen in range(ngens):
        if(ancestral_gen_counter < len(anc_sample_gen) and gen + 1 == anc_sample_gen[ancestral_gen_counter][0]):
            ran_samples = np.random.choice(int(N), int(anc_sample_gen[ancestral_gen_counter][1]), replace=False)
            # while sorting to get diploid chromosomes next to each other isn't strictly necessary,
            # they will be sorted (in reverse order) before simplication anyway, no need to do it here
            ancestral_samples = np.insert(ancestral_samples, len(ancestral_samples), np.concatenate((2*ran_samples + next_id, 2*ran_samples + 1 + next_id)))
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
            simplifier.edges.add_row(0.0, 1.0, p2, next_id+1)
            # Update our dummy variable.
            next_id += 2

        diploids = new_diploids
        new_diploids = diploids + 2*N
        # Let's see if we will do some GC:
        if ((gen+1) % simplifier.gc_interval == 0) or (gen == ngens):
            # If we do GC, we need to reset
            # some variables.  Internally,
            # when msprime simplifies tables,
            # the 2N tips are entries 0 to 2*N-1
            # in the NodeTable, hence the
            # re-assignment of diploids.
            # We can also reset next_id to the
            # length of the current NodeTable,
            # keeping risk of integer overflow
            # to a minimum.
            next_id, diploids, new_diploids, ancestral_samples = simplifier.simplify(gen, ngens, diploids, ancestral_samples)

    return diploids, ancestral_samples
    
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

    ngens = SIMLEN * args.popsize
    simplifier = ARGsimplifier(ngens,args.gc)
    anc_sample_gen = [(ngens * (i + 1) / SIMLEN, max(round(args.popsize / 200), 1)) for i in range(SIMLEN - 2)]
    samples, anc_samples = wf(args.popsize, simplifier, anc_sample_gen, ngens)

    ran_samples = np.random.choice(args.popsize, args.nsam, replace=False)
    nsam_samples = sorted(np.concatenate((2*ran_samples, 2*ran_samples + 1)))
    all_samples = nsam_samples + anc_samples.tolist()
    node_map = msprime.simplify_tables(samples=all_samples, nodes=simplifier.nodes, edges=simplifier.edges)