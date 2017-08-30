import fwdpy11
import fwdpy11.fitness
import fwdpy11.model_params
import numpy as np
import msprime
import time


def evolve_track(rng, pop, params, gc_interval):
    import warnings
    # Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Will throw exception if anything is wrong:
        params.validate()

    from fwdpy11.internal import makeMutationRegions, makeRecombinationRegions
    mm = makeMutationRegions(params.nregions, params.sregions)
    rm = makeRecombinationRegions(params.recregions)

    from .wfarg import evolve_singlepop_regions_track_ancestry, AncestryTracker
    from .argsimplifier import ArgSimplifier
    simplifier = ArgSimplifier(gc_interval)
    atracker = AncestryTracker(pop.N)
    evolve_singlepop_regions_track_ancestry(rng, pop, atracker, simplifier,
                                            params.demography,
                                            params.mutrate_s,
                                            params.recrate, mm, rm,
                                            params.gvalue, params.pself)
    # TODO: need a way to make sure that this
    # is necessary AND to make sure atracker
    # is post-processed!
    simplifier.simplify(atracker)
    return (simplifier,atracker)


def evolve_track_wrapper(popsize=1000, rho=10000.0, mu=1e-2, seed=42,
                         gc_interval=10,
                         dfe=fwdpy11.ConstantS(0, 1, 1, -0.025, 1.0)):
    """
    gc_interval does nothing :)
    """
    if isinstance(dfe, fwdpy11.Sregion) is False:
        raise TypeError("dfe must be a fwdpy11.Sregion")

    if dfe.b != 0.0 or dfe.e != 1.0:
        raise ValueError("DFE beg/end must be 0.0/1.0, repsectively")

    pop = fwdpy11.SlocusPop(popsize)
    recrate = float(rho) / (4.0 * float(popsize))

    pdict = {'rates': (0.0, mu, recrate),
             'nregions': [],
             'sregions': [dfe],
             'recregions': [fwdpy11.Region(0, 1, 1)],
             'gvalue': fwdpy11.fitness.SlocusMult(2.0),
             'demography': np.array([popsize] * 10 * popsize, dtype=np.uint32)
             }

    params = fwdpy11.model_params.SlocusParams(**pdict)
    rng = fwdpy11.GSLrng(seed)
    start_sim = time.time()
    return evolve_track(rng, pop, params, gc_interval)
    # print("got atracker at gen:",atracker.offspring_generation)
    # print(simplifier.nodes)
    # print(len(np.array(atracker.nodes,copy=False)))
    # stop_sim = time.time()

    # Get our nodes and edges from C++ directly

    # start_fudge = time.time()
    # atracker.prep_for_gc();
    # sim_nodes = np.array(atracker.nodes, copy=False)
    # sim_edges = np.array(atracker.edges, copy=False)
    # samples = np.array(atracker.samples, copy=False)
    # # Get sample ids.  Again, better done
    # # on the C++ side
    # # max_gen = sim_nodes['generation'].max()
    # # samples = sim_nodes['id'][np.where(sim_nodes['generation']==max_gen)]
    # # Get the node times, convert to float,
    # # then convert to backwards in time.
    # # This is DUMB and should be handled on the C++
    # # side.
    # # sim_nodes['generation'] -= max_gen
    # # sim_nodes['generation'] *= -1.0
    # stop_fudge = time.time()

    # start_msprime = time.time()
    # n = msprime.NodeTable()
    # flags=np.empty([len(sim_nodes)], dtype=np.uint32)
    # flags.fill(0)
    # is_sample=np.empty([len(samples)], dtype = flags.dtype)
    # is_sample.fill(1)
    # flags[-len(samples):]=is_sample
    # n.set_columns(flags=flags,
    #               # gives type conversion error from uint32 to int32
    #               # without this CAST:
    #               population=sim_nodes['population'], #.astype(np.int32),
    #               time=sim_nodes['generation'])
    # # print(n)
    # e = msprime.EdgesetTable()

    # # This is slow:
    # # for se in sim_edges:
    # #     e.add_row(left=se['left'],
    # #               right=se['right'],
    # #               parent=se['parent'],
    # #               children=(se['child'],))

    # e.set_columns(left=sim_edges['left'],
    #               right=sim_edges['right'],
    #               # CAST
    #               parent=sim_edges['parent'], #.astype(np.int32),
    #               # CAST
    #               children=sim_edges['child'], #.astype(np.int32),
    #               children_length=[1]*len(sim_edges))
    # msprime.sort_tables(nodes=n, edgesets=e)
    # x = msprime.load_tables(nodes=n, edgesets=e)
    # x = x.simplify(samples=samples.tolist())
    # x.dump_tables(nodes=n,edgesets=e)
    # # Based on Peter Ralph's example
    # # at https://github.com/petrelharp/local_pca/blob/master/sims/msp/msp-add-mutation.py
    # msp_rng = msprime.RandomGenerator(seed)
    # mutations = msprime.MutationTable()
    # sites = msprime.SiteTable()
    # mutgen = msprime.MutationGenerator(msp_rng, recrate) # rho = theta
    # mutgen.generate(n, e, sites, mutations)
    # x = msprime.load_tables(nodes=n, edgesets=e, sites=sites, mutations=mutations)
    # stop_msprime = time.time()
    # return {'sim_time': stop_sim - start_sim,
    #         'msprime_time': stop_msprime - start_msprime,
    #         'fudge_time': stop_fudge - start_fudge}
