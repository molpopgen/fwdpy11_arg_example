import msprime
import concurrent.futures
import numpy as np

# You need a function to run


def runsim(args):
    nsam, theta, rho, seed = args
    ts = msprime.simulate(nsam, mutation_rate=theta/4.,
                          recombination_rate=rho/4.,
                          random_seed=seed)
    return ts.tables.sites, ts.tables.mutations, ts.tables.nodes, ts.tables.edges


if __name__ == "__main__":
    np.random.seed(42)
    # Get 100 seeds w/0 replacement from [0,1e6)
    seeds = np.random.choice(range(1000000), 100, replace=False)

    args = [(10, 10., 10., i) for i in seeds]

    with concurrent.futures.ProcessPoolExecutor() as pool:
        futures = {pool.submit(runsim, i) for i in args}
        for fut in concurrent.futures.as_completed(futures):
            sites, mutations, nodes, edges = fut.result()
            # I have msprime 0.6.0.  There is no more
            # load_tables...unsure what to do?
            # The docs for TableCollection don't seem
            # to let me create one.
            # Moreover, this function fails:
            # ts = msprime.load_tables(nodes=nodes,sites=sites,edges=edges,mutations=mutations)
            tc = msprime.TableCollection(1.0)
            tc.nodes.set_columns(time=nodes.time, flags=nodes.flags)
            tc.edges.set_columns(
                parent=edges.parent, child=edges.child, left=edges.left, right=edges.right)
            tc.sites.set_columns(position=sites.position, ancestral_state=sites.ancestral_state,
                                 ancestral_state_offset=sites.ancestral_state_offset)
            tc.mutations.set_columns(site=mutations.site, node=mutations.node,
                                     derived_state=mutations.derived_state, derived_state_offset=mutations.derived_state_offset)
            ts=tc.tree_sequence()
