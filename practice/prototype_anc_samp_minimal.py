# Minimal prototype:
# Does GC at user-specified intervals, ancestral sampling, prior history
import numpy as np
import msprime
import sys
import argparse

# Simulation constants
popsize, ngens, nsam, seed, anc_sample_interval, gc_interval = int(500), int(20*500), int(5), int(42), int(500), int(100)

class ARGsimplifier(object):
    """
    Mimicking the API of our Python
    class to collect simulated
    results and process them via msprime
    """
    def __init__(self):
        self.nodes, self.edges = msprime.NodeTable(), msprime.EdgeTable()
        prior_ts = msprime.simulate(sample_size=2*popsize, Ne=2*popsize, random_seed=seed)
        prior_ts.dump_tables(nodes=self.nodes, edges=self.edges)
        self.nodes.set_columns(flags=self.nodes.flags, population=self.nodes.population, time=self.nodes.time + ngens)

    def simplify(self, generation, samples, anc_samples):
        """
        Details of simplifying.

        :return: length of simplifed node table, which is next_id to use
        """
        diploids = np.arange(2*popsize, dtype=np.uint32)
        if(generation > 0): 
            # Sort and simplify
            msprime.sort_tables(nodes=self.nodes, edges=self.edges)
            # below guards against duplicate node ids when an ancestral sample generation overlaps with a gc generation,
            # ensures that generational samples stay in the same order as in the fw sim and occur *before* ancestral samples in the node table,
            # makes bookkeeping easier during the WF (parents of the next generation are guaranteed to be 0-2N in the node table)
            all_samples = samples.tolist() + [i for i in sorted(anc_samples.tolist()) if i not in samples.tolist()]
            node_map = msprime.simplify_tables(samples=all_samples, nodes=self.nodes, edges=self.edges)
            anc_samples = np.array([node_map[int(node_id)] for node_id in anc_samples])
            # Return length of NodeTable, which can be used as next offspring ID and reset diploids (parent gen to 0, 2N)
        return self.nodes.num_rows, diploids, diploids + self.nodes.num_rows, anc_samples

def wf(simplifier):
    """
    For #popsize diploids, the diploids list contains 2*popsize values.
    For the i-th diploids, diploids[2*i] and diploids[2*i+1]
    represent the two chromosomes.

    We simulate constant popsize here, according to a Wright-Fisher
    scheme:

    1. Pick two parental indexes, each [0,popsize-1]
    2. The two chromosome ids are diploids[2*p] and diploids[2*p+1]
       for each parental index 'p'.
    3. We pass one parental index on to each offspring, according to
       Mendel, which means we swap parental chromosome ids 50% of the time.
    """
    next_id, diploids, new_diploids, ancestral_samples = simplifier.simplify(0, np.empty([0], dtype=np.uint32), np.empty([0], dtype=np.uint32))  #  Based on prior history, this will be the next unique ID to use
    for gen in range(ngens):
        if (((gen + 1) % anc_sample_interval) == 0) and ((gen+1) != ngens): #take 1 ancestral sample every anc_sample_interval generations
            ran_samples = np.random.choice(popsize, 1, replace=False)[0] 
            # while sorting to get diploid chromosomes next to each other isn't strictly necessary,
            # they will be sorted (in reverse order) before simplication anyway, no need to do it here
            ancestral_samples = np.insert(ancestral_samples, len(ancestral_samples), np.array([2*ran_samples + next_id, 2*ran_samples + next_id + 1],dtype=np.uint32))
        # Store nodes for this generation.
        for i in range(2*popsize): simplifier.nodes.add_row(1,(ngens - (gen+1)),0)
        # Pick 2*popsize parents:
        parents = np.random.randint(0, popsize, 2*popsize)
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
        new_diploids = diploids + 2*popsize
        # Let's see if we will do some GC:
        if ((gen+1) % gc_interval == 0) or ((gen+1) == ngens):
            # If we do GC, we need to reset
            # some variables.  Internally,
            # when msprime simplifies tables,
            # the 2*popsize tips are entries 0 to 2*popsize-1
            # in the NodeTable, hence the
            # re-assignment of diploids.
            # We can also reset next_id to the
            # length of the current NodeTable,
            # keeping risk of integer overflow
            # to a minimum.
            next_id, diploids, new_diploids, ancestral_samples = simplifier.simplify(gen, diploids, ancestral_samples)
    return diploids, ancestral_samples
    
if __name__ == "__main__":
    np.random.seed(seed)
    simplifier = ARGsimplifier()
    samples, anc_samples = wf(simplifier)
    ran_samples = np.random.choice(popsize, nsam, replace=False)
    nsam_samples = sorted(np.concatenate((2*ran_samples, 2*ran_samples + 1)))
    all_samples = nsam_samples + anc_samples.tolist()
    node_map = msprime.simplify_tables(samples=all_samples, nodes=simplifier.nodes, edges=simplifier.edges)